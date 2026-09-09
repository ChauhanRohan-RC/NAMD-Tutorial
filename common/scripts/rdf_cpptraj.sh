#!/bin/bash
# =======================================================================
# Bash Wrapper for CPPTRAJ Radial Distribution Function (RDF) g(r)
#
# RUN with tmux:
# $ tmux new-session -d -s name './rdf_cpptraj.sh'
# =======================================================================

# OUTPUT Columns
# 1. r        : radial distance from reference (Å) (bin center)
# 2. g_r      : radial distribution function (normalized probability)
# 3. int_r    : running integral count (coordination number).
#               Avg number of selection_2 atoms within the sphere of radius r

# -------------------------------------------------------------------------------------------------
# HELPER FUNCTION: Find files with a prefix, suffix and an optional number within a given range
# --------------------------------------------------------------------
find_files() {
    local r=() rx="^$2([0-9]+)$3$" n p; for f in "$1"/*; do
        [[ -f "$f" ]] && [[ "${f##*/}" =~ $rx ]] && { n=$((10#${BASH_REMATCH[1]}));
        [[ ( -z "$4" || n -ge $4 ) && ( -z "$5" || n -le $5 ) ]] && {
        [[ "${7:-0}" == 1 ]] && p=$(realpath "$f" 2>/dev/null || readlink -f "$f") || p="$f"; r+=("$p"); }; }
    done; [[ ${#r[@]} -gt 0 ]] && { [[ "${6:-1}" == 1 ]] && printf "%s\n" "${r[@]}" | sort -V || printf "%s\n" "${r[@]}"; }
}


# ------------------------------------
# INPUT Params (Search TODO)
# ------------------------------------

module load ambertools

# Execution Settings
CPPTRAJ_BIN="cpptraj.OMP"          # Executable name (often cpptraj.OMP for OpenMP builds)
OMP_THREADS=16                     # [only with cpptraj.OMP] Number of CPU threads for OpenMP builds, 1 to disable
DEBUG_MODE=1                       # 1 = ON (Keep raw cpptraj outputs and logs in /tmp), 0 = OFF (Delete them)

# RAM Cache Settings (Optional Feature)
USE_RAM_CACHE=1                # 1 = ON (Process sequentially in RAM), 0 = OFF (Process all from disk simultaneously)
RAM_CACHE_DIR="/tmp"           # Target directory for temporary storage (tmpfs/RAM disk)

# File Inputs
PSF_FILE="../../common/amyl_wb.psf"   # TODO: input structure (.psf) file
mapfile -t DCD_FILES < <(find_files ".." "amyl_wb_eq" ".dcd" "" "" 1 0)

# Selections (AMBER Mask Syntax)
SOLUTE_MASK="!(:TIP3,SOD,CLA)&!(@H*)"           # TODO: Reference subsystem (usually heavy protein atoms)
SOLVENT_MASK=":TIP3@OH2"                        # TODO: distribution subsystem

# RDF Parameters
MAX_RADIUS=10.0                # Maximum distance for RDF calculation (Å)
BIN_SIZE=0.05                  # Distance bin size (Å)
USE_PBC=1                      # 1 = ON (Default), 0 = OFF (adds 'noimage' flag)
USE_VOLUME=1                   # 1 = ON (NPT density normalization), 0 = OFF
FRAME_STEP=1                   # Read every Nth frame

# Output
OUT_FILE="rdf.csv"             # TODO: output file
OUT_DELIMITER=" "
COMMENT_TOKEN="#"


# =======================================================================
# SIGNAL TRAP FOR CLEANUP
# =======================================================================
# This ensures temporary files are safely deleted if you Ctrl+C or kill the script
cleanup_on_interrupt() {
    echo -e "\n[ WARN   ] Script interrupted (SIGINT/SIGTERM)! Cleaning up temporary files..."
    
    # Safely delete script and active dcd from memory
    [[ -n "$TMP_SCRIPT" && -f "$TMP_SCRIPT" ]] && rm -f "$TMP_SCRIPT"
    [[ -n "$run_dcd" && -f "$run_dcd" ]] && rm -f "$run_dcd"
    
    if [ "${DEBUG_MODE:-0}" -eq 1 ]; then
        echo "[ DEBUG  ] Interrupt caught. Keeping raw logs and subset CSVs in /tmp for inspection."
    else
        [[ -n "$TMP_LOG" && -f "$TMP_LOG" ]] && rm -f "$TMP_LOG"
        
        # Safely loop and delete output array elements
        if [[ ${#TMP_OUT_LIST[@]} -gt 0 ]]; then
            for tmp_out in "${TMP_OUT_LIST[@]}"; do
                [[ -n "$tmp_out" && -f "$tmp_out" ]] && rm -f "$tmp_out"
            done
        fi
    fi
    exit 1
}
trap cleanup_on_interrupt SIGINT SIGTERM


# =======================================================================
# MAIN
# =======================================================================

# Force standard C locale to prevent awk from truncating floats (e.g., expecting commas instead of dots)
export LC_ALL=C

# Benchmark Trackers
SCRIPT_START_TIME=$SECONDS
TOTAL_IO_TIME=0
TOTAL_COMPUTE_TIME=0
TOTAL_AVERAGING_TIME=0
TOTAL_POST_TIME=0

echo "==================================================================="
echo "[  INIT  ] Starting CPPTRAJ Dynamic RDF Calculator"
echo "[ CONFIG ] Structure     : $PSF_FILE"
echo "[ CONFIG ] Solute Mask   : $SOLUTE_MASK"
echo "[ CONFIG ] Solvent Mask  : $SOLVENT_MASK"
echo "[ CONFIG ] Radius (Max)  : $MAX_RADIUS A | Bin Size: $BIN_SIZE A"
echo "[ CONFIG ] Frame Step    : $FRAME_STEP"
echo "[ CONFIG ] Periodic Bound: $( [ "$USE_PBC" -eq 1 ] && echo "Enabled" || echo "Disabled" )"
echo "[ CONFIG ] RAM Cache     : $( [ "$USE_RAM_CACHE" -eq 1 ] && echo "Enabled ($RAM_CACHE_DIR)" || echo "Disabled" )"
echo "[ CONFIG ] Debug Mode    : $( [ "$DEBUG_MODE" -eq 1 ] && echo "Enabled (Logs kept)" || echo "Disabled" )"
echo "==================================================================="

# Check if CPPTRAJ is in PATH
if ! command -v "$CPPTRAJ_BIN" &> /dev/null; then
    echo "[ ERROR  ] Executable '$CPPTRAJ_BIN' could not be found."
    exit 1
fi

# Configure OpenMP Threading
if [ "$OMP_THREADS" -gt 1 ]; then
    export OMP_NUM_THREADS=$OMP_THREADS
    export OMP_PROC_BIND=false
    export OMP_PLACES=cores
    export OMP_MAX_ACTIVE_LEVELS=1
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1

    if ! echo "quit" | "$CPPTRAJ_BIN" 2>&1 | grep -qi "OpenMP"; then
        echo "[ WARN   ] OMP_THREADS > 1, but '$CPPTRAJ_BIN' does not appear to be an OpenMP build!"
    else
        echo "[ SYSTEM ] OpenMP threading enabled using $OMP_THREADS threads."
    fi
else
    export OMP_NUM_THREADS=1
fi


# Format Optional Flags
VOL_FLAG=""
PBC_FLAG=""
[ "$USE_VOLUME" -eq 1 ] && VOL_FLAG="volume"
[ "$USE_PBC" -eq 0 ] && PBC_FLAG="noimage"

cpptraj_exit_code=0


# =======================================================================
# EXECUTION & DYNAMIC SCRIPT GENERATION
# =======================================================================

if [ "$USE_RAM_CACHE" -eq 1 ]; then
    echo "-------------------------------------------------------------------"
    echo "[ SYSTEM ] RAM Cache enabled. Processing trajectories sequentially."
    echo "-------------------------------------------------------------------"

    TMP_OUT_LIST=()
    FRAME_COUNTS=()
    WEIGHTED_AVERAGE=1

    for i in "${!DCD_FILES[@]}"; do
        dcd="${DCD_FILES[$i]}"
        dcd_base=$(basename "$dcd")
        
        # Use mktemp to ensure unique filenames without collisions
        tmp_out=$(mktemp -t tmp_rdf_out_XXXXXX.dat)
        TMP_OUT_LIST+=("$tmp_out")

        if [ ! -f "$dcd" ]; then
            echo "[ ERROR  ] Trajectory file not found: $dcd"
            exit 1
        fi

        # Safely convert to KB to avoid 32-bit bash integer overflow for files > 2GB
        file_size_kb=$(( $(stat -c%s "$dcd" 2>/dev/null || stat -f%z "$dcd" 2>/dev/null) / 1024 ))
        free_space_kb=$(df -k "$RAM_CACHE_DIR" 2>/dev/null | awk 'NR==2 {print $4}')
        buffer_kb=51200 # 50 MB safety buffer

        cached=0
        if [ -n "$file_size_kb" ] && [ -n "$free_space_kb" ] && [ "$free_space_kb" -gt $((file_size_kb + buffer_kb)) ]; then
            run_dcd=$(mktemp "$RAM_CACHE_DIR/tmp_dcd_XXXXXX.dcd")
            echo "[ INFO   ] [$((i+1))/${#DCD_FILES[@]}] Copying $dcd_base to RAM ($run_dcd) ..."
            
            io_start=$SECONDS
            cp "$dcd" "$run_dcd"
            TOTAL_IO_TIME=$((TOTAL_IO_TIME + SECONDS - io_start))
            
            cached=1
        else
            echo "[ WARN   ] [$((i+1))/${#DCD_FILES[@]}] Insufficient space in RAM. Falling back to disk: $dcd_base"
            run_dcd="$dcd"
        fi

        # Generate CPPTRAJ input script for the single trajectory
        TMP_SCRIPT=$(mktemp -t cpptraj_rdf_XXXXXX.in)
        TMP_LOG=$(mktemp -t cpptraj_log_XXXXXX.out)

        cat << EOF > "$TMP_SCRIPT"
parm $PSF_FILE
trajin $run_dcd 1 last $FRAME_STEP
radial $tmp_out $BIN_SIZE $MAX_RADIUS $SOLVENT_MASK $SOLUTE_MASK intrdf $tmp_out $PBC_FLAG $VOL_FLAG
go
quit
EOF

        echo "[ EXEC   ] Processing $dcd_base with $CPPTRAJ_BIN ..."
        
        comp_start=$SECONDS
        "$CPPTRAJ_BIN" -i "$TMP_SCRIPT" > "$TMP_LOG" 2>&1
        local_exit=$?
        TOTAL_COMPUTE_TIME=$((TOTAL_COMPUTE_TIME + SECONDS - comp_start))
        
        # Bulletproof extraction: Isolate 'processed X frames' then extract only the numbers
        processed_frames=$(grep -oE 'processed [0-9]+ frames' "$TMP_LOG" | grep -oE '[0-9]+' | tail -n 1)

        # Fallback in case of formatting mismatch or missing logs
        if ! [[ "$processed_frames" =~ ^[0-9]+$ ]] || [ "$processed_frames" -eq 0 ]; then
            echo "[ WARN   ] Strict frame extraction failed for $dcd_base."
            echo "[ WARN   ] Block averaging will completely fallback to non-weighted arithmetic mean."
            WEIGHTED_AVERAGE=0
        else
            FRAME_COUNTS+=("$processed_frames")
            echo "[ INFO   ] Success: Computed $processed_frames frames from $dcd_base."
        fi

        if [ $local_exit -ne 0 ]; then
            echo "[ ERROR  ] CPPTRAJ encountered an error on $dcd_base (exit: $local_exit). Check raw logs."
            cpptraj_exit_code=$local_exit
        fi

        # Cleanup RAM Cache & Temp Files
        rm -f "$TMP_SCRIPT"
        
        if [ "$DEBUG_MODE" -eq 1 ]; then
            echo "[ DEBUG  ] Kept run log: $TMP_LOG"
        else
            rm -f "$TMP_LOG"
        fi

        if [ "$cached" -eq 1 ]; then
            echo "[ INFO   ] Deleting temporary DCD from RAM ..."
            rm -f "$run_dcd"
        fi
        echo "-------------------------------------------------------------------"
    done

    # Block average all temporary subset outputs into the final output
    if [ $cpptraj_exit_code -eq 0 ]; then
        avg_start=$SECONDS
        
        if [ "$WEIGHTED_AVERAGE" -eq 1 ]; then
            echo "[ SYSTEM ] Performing mathematical FRAME-WEIGHTED block averaging..."
            # Safely pass array to awk as string space separated
            awk -v weights="${FRAME_COUNTS[*]}" '
            BEGIN {
                split(weights, w, " ")
                total_frames = 0
                for (i in w) {
                    total_frames += (w[i] + 0)
                }
                file_idx = 0
            }
            FNR == 1 {
                file_idx++
            }
            /^#/ { 
                next 
            }
            # Ensure line has at least 3 columns of data to avoid math errors on blank lines
            NF >= 3 {
                r[FNR] = $1
                weight = (w[file_idx] != "") ? (w[file_idx] + 0) : 1
                g[FNR] += ($2 + 0.0) * weight
                int_r[FNR] += ($3 + 0.0) * weight
                if (FNR > max_fnr) max_fnr = FNR
            }
            END {
                for (i=1; i<=max_fnr; i++) {
                    if (r[i] != "") {
                        if (total_frames > 0) {
                            printf "%s %.6f %.6f\n", r[i], g[i]/total_frames, int_r[i]/total_frames
                        } else {
                            printf "%s 0.000000 0.000000\n", r[i]
                        }
                    }
                }
            }' "${TMP_OUT_LIST[@]}" > "$OUT_FILE"
            echo "[ INFO   ] Weighted block averaging complete. Total frames processed: $(IFS=+; echo "$((${FRAME_COUNTS[*]}))")"

        else
            echo "[ SYSTEM ] Performing standard UNWEIGHTED block averaging (Fallback Mode)..."
            awk '
            /^#/ { 
                next 
            }
            NF >= 3 {
                r[FNR] = $1
                g[FNR] += ($2 + 0.0)
                int_r[FNR] += ($3 + 0.0)
                count[FNR]++
                if (FNR > max_fnr) max_fnr = FNR
            }
            END {
                for (i=1; i<=max_fnr; i++) {
                    if (count[i] > 0 && r[i] != "") {
                        printf "%s %.6f %.6f\n", r[i], g[i]/count[i], int_r[i]/count[i]
                    }
                }
            }' "${TMP_OUT_LIST[@]}" > "$OUT_FILE"
            echo "[ INFO   ] Unweighted block averaging complete."
        fi
        
        TOTAL_AVERAGING_TIME=$((SECONDS - avg_start))
    fi
    
    # Cleanup temp outputs
    if [ "$DEBUG_MODE" -eq 1 ]; then
        echo "[ DEBUG  ] Kept raw output CSVs: ${TMP_OUT_LIST[*]}"
    else
        rm -f "${TMP_OUT_LIST[@]}"
    fi

else
    # ORIGINAL DISK-BASED BATCH EXECUTION
    echo "-------------------------------------------------------------------"
    echo "[ SYSTEM ] Processing all trajectories in a single batch (Disk Mode)."
    echo "-------------------------------------------------------------------"
    TMP_SCRIPT=$(mktemp -t cpptraj_rdf_XXXXXX.in)
    
    cat << EOF > "$TMP_SCRIPT"
parm $PSF_FILE
EOF

    for dcd in "${DCD_FILES[@]}"; do
        if [ ! -f "$dcd" ]; then
            echo "[ ERROR  ] Trajectory file not found: $dcd"
            rm -f "$TMP_SCRIPT"
            exit 1
        fi
        echo "trajin $dcd 1 last $FRAME_STEP" >> "$TMP_SCRIPT"
    done

    cat << EOF >> "$TMP_SCRIPT"
radial $OUT_FILE $BIN_SIZE $MAX_RADIUS $SOLVENT_MASK $SOLUTE_MASK intrdf $OUT_FILE $PBC_FLAG $VOL_FLAG
go
quit
EOF

    echo "[ EXEC   ] Launching $CPPTRAJ_BIN engine..."
    
    comp_start=$SECONDS
    "$CPPTRAJ_BIN" -i "$TMP_SCRIPT"
    cpptraj_exit_code=$?
    TOTAL_COMPUTE_TIME=$((SECONDS - comp_start))
    
    rm -f "$TMP_SCRIPT"
fi


# =======================================================================
# POST-PROCESSING
# =======================================================================

if [ -f "$OUT_FILE" ] && [ $cpptraj_exit_code -eq 0 ]; then
    post_start=$SECONDS
    TMP_OUT=$(mktemp)

    cat << EOF > "$TMP_OUT"
${COMMENT_TOKEN} ===================================================================
${COMMENT_TOKEN} CPPTRAJ Radial Distribution Function (RDF) & Integral
${COMMENT_TOKEN} ===================================================================
${COMMENT_TOKEN} Structure File : $PSF_FILE
${COMMENT_TOKEN} Trajectories   : ${DCD_FILES[*]}
${COMMENT_TOKEN} -----------------------------------------------
${COMMENT_TOKEN} Solute Mask    : $SOLUTE_MASK    (Selection-1, reference subsystem)
${COMMENT_TOKEN} Solvent Mask   : $SOLVENT_MASK   (Selection-2, distributed around reference)
${COMMENT_TOKEN} -----------------------------------------------
${COMMENT_TOKEN} Max Radius     : $MAX_RADIUS Å
${COMMENT_TOKEN} Bin Size       : $BIN_SIZE Å
${COMMENT_TOKEN} Periodic Bound : $( [ "$USE_PBC" -eq 1 ] && echo "Enabled" || echo "Disabled (noimage)" )
${COMMENT_TOKEN} Volume Norm    : $( [ "$USE_VOLUME" -eq 1 ] && echo "Enabled" || echo "Disabled" )
${COMMENT_TOKEN} Frame Step     : $FRAME_STEP
${COMMENT_TOKEN} RAM Cache Mode : $( [ "$USE_RAM_CACHE" -eq 1 ] && echo "Enabled" || echo "Disabled" )
${COMMENT_TOKEN} -----------------------------------------------
${COMMENT_TOKEN} OUTPUT
${COMMENT_TOKEN} => r        : radial distance from reference (Å) (bin center)
${COMMENT_TOKEN} => g_r      : radial distribution function RDF (normalized probability)
${COMMENT_TOKEN} => int_r    : running integral count (coordination number). Avg number of selection-2 atoms within the sphere of radius r
${COMMENT_TOKEN} -----------------------------------------------
EOF

    out_header="r${OUT_DELIMITER}g_r${OUT_DELIMITER}int_r"

    awk -v d="$OUT_DELIMITER" -v c="$COMMENT_TOKEN" -v h="$out_header" '
    BEGIN {
        print h
    }
    /^#/ {
        print c " [ARCHIVE] Original CPPTRAJ Header: " $0
        print c "============================================================"
        next
    }
    {
        # In batch mode CPPTRAJ output, force 6 decimals here too
        printf "%s%s%.6f%s%.6f\n", $1, d, $2, d, $3
    }' "$OUT_FILE" >> "$TMP_OUT"

    mv "$TMP_OUT" "$OUT_FILE"
    
    TOTAL_POST_TIME=$((SECONDS - post_start))
fi


# =======================================================================
# FINALIZING & BENCHMARKING
# =======================================================================
if [ $cpptraj_exit_code -eq 0 ] && [ -f "$OUT_FILE" ]; then
    echo "-------------------------------------------------------------------"
    echo "[ SUCCESS ]"
else
    echo "-------------------------------------------------------------------"
    echo "[ ERROR   ] CPPTRAJ encountered an error (exit_code $cpptraj_exit_code) OR output file was not generated."
fi

if [ -f "$OUT_FILE" ]; then
    echo "[ OUTPUT ] Output File: \"$OUT_FILE\""
fi

TOTAL_WALL_TIME=$((SECONDS - SCRIPT_START_TIME))

echo "==================================================================="
echo "[ BENCHMARK ] Execution Time Summary"
echo "-------------------------------------------------------------------"
if [ "$USE_RAM_CACHE" -eq 1 ]; then
    printf " %-26s | %5d seconds\n" "Total I/O (RAM Cache)" "$TOTAL_IO_TIME"
    printf " %-26s | %5d seconds\n" "Total CPPTRAJ Compute" "$TOTAL_COMPUTE_TIME"
    printf " %-26s | %5d seconds\n" "Block Averaging" "$TOTAL_AVERAGING_TIME"
else
    printf " %-26s | %5d seconds\n" "Total CPPTRAJ Compute" "$TOTAL_COMPUTE_TIME"
fi
printf " %-26s | %5d seconds\n" "Post-Processing" "$TOTAL_POST_TIME"
echo "-------------------------------------------------------------------"
printf " %-26s | %5d seconds\n" "Total Overall Wall Time" "$TOTAL_WALL_TIME"
echo "==================================================================="
