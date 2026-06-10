#!/bin/bash

# --------------------------------------------------------------------
# HELPER FUNCTION: Find files with a prefix, suffix and an optional number within a given range
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
find_files() {
    local dir_path="$1"
    local prefix="$2"
    local suffix="$3"
    local min_num="${4:-}"
    local max_num="${5:-}"
    local sort_natural="${6:-1}"
    local return_abs_path="${7:-0}"

    # Bash regex (no quotes around the regex string when using in [[ =~ ]])
    local regex="^${prefix}([0-9]+)${suffix}$"
    local result_list=()

    # Iterate through directory safely handling spaces in filenames
    for f in "$dir_path"/*; do
        # Skip if not a regular file (also handles empty directory globbing)
        [[ -f "$f" ]] || continue

        local filename
        filename=$(basename -- "$f")

        # Check regex match
        if [[ $filename =~ $regex ]]; then
            # Extract the capture group (the number)
            local num="${BASH_REMATCH[1]}"

            # CRITICAL: Force base-10 evaluation.
            # Otherwise, Bash treats numbers with leading zeros (e.g., "08") as invalid octal numbers.
            local num_dec=$((10#$num))

            # Range checks
            if [[ ( -z "$min_num" || $num_dec -ge $min_num ) && ( -z "$max_num" || $num_dec -le $max_num ) ]]; then

                local fpath="$f"
                if [[ "$return_abs_path" -eq 1 ]]; then
                    # Get absolute path (realpath is standard on most modern Linux distros)
                    fpath=$(realpath "$f" 2>/dev/null || readlink -f "$f")
                fi

                result_list+=("$fpath")
            fi
        fi
    done

    # Output the results, optionally passing through a natural version sort
    if [[ ${#result_list[@]} -gt 0 ]]; then
        if [[ "$sort_natural" -eq 1 ]]; then
            printf "%s\n" "${result_list[@]}" | sort -V
        else
            printf "%s\n" "${result_list[@]}"
        fi
    fi
}