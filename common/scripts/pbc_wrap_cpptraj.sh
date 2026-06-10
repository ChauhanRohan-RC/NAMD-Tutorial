#!/bin/bash

# ============================================================================
# Wrap trajectories using cpptraj (AMBERTOOLS) about Periodic Boundary COndition (PBC)
# -> Optimal for large DCD files. For single frames, use pbc_wrap.py
# ---------------------------------------------------------------------------
# => Extremely efficient and fast
# => Optimal for big DCD files
# => DOES NOT Support writing to binary .coor and .restart.coor.
#    Use pbc_wrap.py (MDAnalysis) for those purposes.
# ---------------------------------------------------------------------------
# USGAE: search for TODO
# ============================================================================

# -------------------------------------------------------------------------------------------------
# HELPER FUNCTION: Find files with a prefix, suffix and an optional number within a given range
# USAGE: mapfile -t input_traj_files < <(find_files "." "a" ".dcd" "" "" 1 0)
#
# Arguments:
#   1: dir_path        : directory path to search
#   2: prefix          : file name prefix
#   3: suffix          : file name suffix
#   4: min_num         : optional min range (inclusive). Leave empty ("") for none.
#   5: max_num         : optional max range (inclusive). Leave empty ("") for none.
#   6: sort_natural    : 1 for natural sort (default), 0 for no sort
#   7: return_abs_path : 1 for absolute path, 0 for relative path (default)
# --------------------------------------------------------------------
find_files() {
    local r=() rx="^$2([0-9]+)$3$" n p; for f in "$1"/*; do
        [[ -f "$f" ]] && [[ "${f##*/}" =~ $rx ]] && { n=$((10#${BASH_REMATCH[1]}));
        [[ ( -z "$4" || n -ge $4 ) && ( -z "$5" || n -le $5 ) ]] && {
        [[ "${7:-0}" == 1 ]] && p=$(realpath "$f" 2>/dev/null || readlink -f "$f") || p="$f"; r+=("$p"); }; }
    done; [[ ${#r[@]} -gt 0 ]] && { [[ "${6:-1}" == 1 ]] && printf "%s\n" "${r[@]}" | sort -V || printf "%s\n" "${r[@]}"; }
}

# ------------------------------------------------------------------------------------------------------



input_psf_file="../common/amyld_wb.psf"   # TODO: input structure (.psf) file

# TODO: input array of DCD trajectory files
#input_traj_files=(
#    "amyld_wb_eq50.dcd"
#    "amyld_wb_eq51.dcd"
#)
mapfile -t input_traj_files < <(find_files "." "amyld_wb_eq" ".dcd" "" "" 1 0)

input_xsc_file=""  # [OPTIONAL] only used for non-dcd files lacking box information

overwrite_input_files="false"  # overwrite input traj file(s) with their wrapped version on success
out_dir=""               # [OPTIONAL] output folder. ONLY WHEN overwrite is false. Defaults to the same folder as each of the input file
out_prefix=""            # [OPTIONAL] output file prefix (added to input file name)
out_suffix=".wrapped"    # TODO: output file suffix (added to input file name)

# ---- Wrapping Algorithm ------
# use cpptraj's highly optimized autoimage (Recommended).
# Set to "false" to manually define the centering and wrapping steps below.
use_autoimage="true"    # [true/false]

# ----- Manual Wrapping Settings (Ignored if use_autoimage="true") -----
# CRITICAL: Must use AMBER mask syntax!
# Examples: "*@CA" (all C-alphas), ":1-100" (residues 1 to 100), ":PRO" (all proteins)
center_atom_selection=":PRO"
center_method="mass"    # Centering method: "mass" (weighted) or "" (leave empty for geometry)
wrap_compound="byres"   # Wrapping compound: "byres" (keep residues intact) or "bymol" (keep molecules intact)

# Execution options
use_cuda="false"         # uSe cpptraj.cuda if available
verbose="true"          # show cpptraj output


# ==============================================================================
# MAIN
# ==============================================================================



BOX_X="" ; BOX_Y="" ; BOX_Z=""
BOX_ALPHA="" ; BOX_BETA="" ; BOX_GAMMA=""

if [[ -n "$input_xsc_file" && -f "$input_xsc_file" ]]; then
    echo "Parsing box dimensions from $input_xsc_file..."
    # Extract the last line that doesn't start with '#'
    last_line=$(grep -v "^#" "$input_xsc_file" | tail -n 1)

    # Pass the line to awk to do the heavy math
    BOX_INFO=$(echo "$last_line" | awk '{
        # Define an arccosine function (awk only has atan2 natively)
        function acos(x) {
            if (x <= -1) return 3.141592653589793;
            if (x >= 1) return 0.0;
            return atan2(sqrt(1-x*x), x);
        }

        # NAMD XSC format vectors
        ax=$2; ay=$3; az=$4;
        bx=$5; by=$6; bz=$7;
        cx=$8; cy=$9; cz=$10;

        # Lengths (Pythagorean theorem)
        a = sqrt(ax*ax + ay*ay + az*az);
        b = sqrt(bx*bx + by*by + bz*bz);
        c = sqrt(cx*cx + cy*cy + cz*cz);

        # Angles (Dot product, converted from radians to degrees)
        rad2deg = 180.0 / 3.141592653589793;
        alpha = acos((bx*cx + by*cy + bz*cz) / (b * c)) * rad2deg;
        beta  = acos((ax*cx + ay*cy + az*cz) / (a * c)) * rad2deg;
        gamma = acos((ax*bx + ay*by + az*bz) / (a * b)) * rad2deg;

        # Print formatted output for Bash to read
        printf "%.4f %.4f %.4f %.4f %.4f %.4f", a, b, c, alpha, beta, gamma
    }')

    # Read the awk output back into Bash variables
    read BOX_X BOX_Y BOX_Z BOX_ALPHA BOX_BETA BOX_GAMMA <<< "$BOX_INFO"
    echo "Extracted -> L: ($BOX_X, $BOX_Y, $BOX_Z) Angles: ($BOX_ALPHA, $BOX_BETA, $BOX_GAMMA)"
fi


# -----------------------------
# MAIN EXECUTION LOOP
# -----------------------------
echo "------------------------------------------------"
echo " => Starting cpptraj Batch Wrapping Pipeline "

cmd="cpptraj"
TMP_SCRIPT="tmp_cpptraj_wrap.in"

if [[ "$use_cuda" == "true" ]] && type -P "${cmd}.cuda" >/dev/null 2>&1; then
    cmd="${cmd}.cuda"
    echo " => Using CUDA binary: ${cmd}"
else
    echo " => Using binary: ${cmd}"
fi
echo "------------------------------------------------"

for in_file in "${input_traj_files[@]}"; do
    # 1. Precondition Check
    if [[ ! -f "$in_file" ]]; then
        echo "WARNING: File '$in_file' not found. Skipping."
        continue
    fi

    # 2. Generate Output Filename
    dirname=$(dirname -- "$in_file")
    filename=$(basename -- "$in_file")
    extension="${filename##*.}"
    stem="${filename%.*}"

    if [[ "$overwrite_input_files" == "true" && -z "$out_suffix" && -z "$out_prefix" ]]; then
        out_suffix=".wrapped"   # default suffix when replacing files
        echo "INFO: Using default OUTPUT_SUFFIX = ${out_suffix}"
    fi

    # Determine target directory
    if [[ "$overwrite_input_files" != "true" && -n "$out_dir" ]]; then
        target_dir="$out_dir"
        mkdir -p "$target_dir"
    else
        target_dir="$dirname"
    fi

    out_file="${target_dir}/${out_prefix}${stem}${out_suffix}.${extension}"
    # Check if out_file == file
    if [[ "$in_file" -ef "$out_file" ]]; then
        echo "ERROR: OUTPUT and INPUT file paths are the same -> '${in_file}'"
        continue
    fi

    if [[ -f "$out_file" ]]; then
        echo "ERROR: output file already exists -> '${out_file}'"
        continue
    fi

    echo -e "\n[+] Processing: '${in_file}' -> '${out_file}'"

    echo "parm $input_psf_file" > "$TMP_SCRIPT"
    echo "trajin $in_file" >> "$TMP_SCRIPT"

    # If it is NOT a DCD file, explicitly set the box dimensions using the parsed variables
    # ${extension,,} converts the extension to lowercase for safe comparison
    if [[ "${extension,,}" != "dcd" ]]; then
        if [[ -z "$BOX_X" ]]; then
            echo "CRITICAL ERROR: Non-DCD file detected ($in_file) but no valid XSC file parsed. Exiting."
            rm -f "$TMP_SCRIPT"
            exit 1
        fi
        echo "box x $BOX_X y $BOX_Y z $BOX_Z alpha $BOX_ALPHA beta $BOX_BETA gamma $BOX_GAMMA" >> "$TMP_SCRIPT"
    fi

    # 4. Apply Wrapping Logic
    if [[ "$use_autoimage" == "true" ]]; then
        echo "autoimage" >> "$TMP_SCRIPT"
    else
        # Manual Centering and Imaging
        echo "center $center_atom_selection $center_method origin" >> "$TMP_SCRIPT"
        # 'familiar' applies the standard truncated/orthogonal box geometry
        echo "image origin center familiar $wrap_compound" >> "$TMP_SCRIPT"
    fi

    # 5. Output and Run
    echo "trajout $out_file" >> "$TMP_SCRIPT"
    echo "run" >> "$TMP_SCRIPT"
    echo "quit" >> "$TMP_SCRIPT"

    # Execute cpptraj
    if [[ "$verbose" == "true" ]]; then
        "$cmd" -i "$TMP_SCRIPT"
    else
        "$cmd" -i "$TMP_SCRIPT" > /dev/null 2>&1
    fi

    # Check if the command succeeded
    if [[ $? -eq 0 ]]; then
        echo "SUCCESS: Wrapped trajectory saved."

        # Overwriting input trajectory file
        if [[ "$overwrite_input_files" == "true" ]]; then
            echo "------------------------------------------------------------"
            echo "=> OVER-WRITING input trajectory file with wrapped version"
            echo "   '${out_file}' -> '${in_file}'"

            mv "${out_file}" "${in_file}"
            echo "------------------------------------------------------------"
        fi
    else
        echo "FAILED: cpptraj encountered an error. Run manually without > /dev/null to debug."
    fi

    # Clean up the temporary script for the next iteration
    rm -f "$TMP_SCRIPT"

done

echo "Wrapping Finished"





# --------------------------------------------------------------------
# FULL IMPLEMENTATION: Find files with a prefix, suffix and an optional number within a given range
#
# USAGE EXAMPLE:
#-----------------
# mapfile -t files_arr < <(find_files "." "a" ".dcd" "" "" 1 0)
# echo "Found ${#files_arr[@]} files to process."
# for f in "${files_arr[@]}"; do
#     echo " -> $f"
# done
#
# Arguments:
#   1: dir_path        : directory path to search
#   2: prefix          : file name prefix
#   3: suffix          : file name suffix
#   4: min_num         : optional min range (inclusive). Leave empty ("") for none.
#   5: max_num         : optional max range (inclusive). Leave empty ("") for none.
#   6: sort_natural    : 1 for natural sort (default), 0 for no sort
#   7: return_abs_path : 1 for absolute path, 0 for relative path (default)
# --------------------------------------------------------------------

#find_files() {
#    local dir_path="$1"
#    local prefix="$2"
#    local suffix="$3"
#    local min_num="${4:-}"
#    local max_num="${5:-}"
#    local sort_natural="${6:-1}"
#    local return_abs_path="${7:-0}"
#
#    # Bash regex (no quotes around the regex string when using in [[ =~ ]])
#    local regex="^${prefix}([0-9]+)${suffix}$"
#    local result_list=()
#
#    # Iterate through directory safely handling spaces in filenames
#    for f in "$dir_path"/*; do
#        # Skip if not a regular file (also handles empty directory globbing)
#        [[ -f "$f" ]] || continue
#
#        local filename
#        filename=$(basename -- "$f")
#
#        # Check regex match
#        if [[ $filename =~ $regex ]]; then
#            # Extract the capture group (the number)
#            local num="${BASH_REMATCH[1]}"
#
#            # CRITICAL: Force base-10 evaluation.
#            # Otherwise, Bash treats numbers with leading zeros (e.g., "08") as invalid octal numbers.
#            local num_dec=$((10#$num))
#
#            # Range checks
#            if [[ ( -z "$min_num" || $num_dec -ge $min_num ) && ( -z "$max_num" || $num_dec -le $max_num ) ]]; then
#
#                local fpath="$f"
#                if [[ "$return_abs_path" -eq 1 ]]; then
#                    # Get absolute path (realpath is standard on most modern Linux distros)
#                    fpath=$(realpath "$f" 2>/dev/null || readlink -f "$f")
#                fi
#
#                result_list+=("$fpath")
#            fi
#        fi
#    done
#
#    # Output the results, optionally passing through a natural version sort
#    if [[ ${#result_list[@]} -gt 0 ]]; then
#        if [[ "$sort_natural" -eq 1 ]]; then
#            printf "%s\n" "${result_list[@]}" | sort -V
#        else
#            printf "%s\n" "${result_list[@]}"
#        fi
#    fi
#}