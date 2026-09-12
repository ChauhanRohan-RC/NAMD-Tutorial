#!/usr/bin/env python3

import os
import glob

#----------------------------------------------------------
# Simple script to drop all zero-columns from data files
#----------------------------------------------------------

INPUT_FILES: list[str] = glob.glob('*.energy.csv')     # list of input data files

ZERO_THRESHOLD: float = 0.0     # a value is considered 0 if abs(value) <= threshold.
                                  # set exactly to 0.0 to only drop exact zeroes

OUTPUT_SUFFIX = '.filtered.csv'
COMMENT_TOKEN = '#'
OUTPUT_DELIMITER = ' '         



# ------------------------------
# MAIN
#-------------------------------

# Terminal colors for pretty logging
class LogStyle:
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    RESET = '\033[0m'

def log_info(msg):
    print(f"{LogStyle.CYAN}{LogStyle.BOLD}[INFO]{LogStyle.RESET} {msg}")

def log_success(msg):
    print(f"{LogStyle.GREEN}{LogStyle.BOLD}[SUCCESS]{LogStyle.RESET} {msg}")

def log_warn(msg):
    print(f"{LogStyle.YELLOW}{LogStyle.BOLD}[SKIP]{LogStyle.RESET} {msg}")

def process_csv(filepath):
    comments = []
    data = []

    # 1. Read file and separate comments from data
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith(COMMENT_TOKEN):
                comments.append(line)
            else:
                line_stripped = line.strip()
                if line_stripped:
                    # Splits by arbitrary whitespace if no delimiter is provided to split()
                    data.append(line_stripped.split())

    if not data:
        log_warn(f"'{filepath}': No data rows found.")
        return

    # 2. Identify which columns to keep
    # Get the max number of columns in case rows are jagged
    max_cols = max(len(row) for row in data)
    cols_to_keep = []

    for col_idx in range(max_cols):
        has_numeric = False
        has_non_zero = False

        for row in data:
            if col_idx < len(row):
                val = row[col_idx]
                try:
                    f_val = float(val)
                    has_numeric = True
                    # If any numeric value is not mathematically zero, we keep the column
                    if abs(f_val) > ZERO_THRESHOLD:
                        has_non_zero = True
                except ValueError:
                    # Ignore headers or string data. We don't fail here.
                    pass

        # KEEP logic:
        # - If it has non-zero numbers, keep it.
        # - If it has NO numbers at all (purely string categories/text), keep it.
        # Drop ONLY if it has numbers and ALL numbers evaluate to 0.0
        if has_non_zero or not has_numeric:
            cols_to_keep.append(col_idx)

    # 3. Write out to new file
    base_name, ext = os.path.splitext(filepath)
    out_filepath = f"{base_name}{OUTPUT_SUFFIX}"

    with open(out_filepath, 'w') as f:
        # Write preserved comments
        for comment in comments:
            f.write(comment)

        # Write filtered data rows
        for row in data:
            filtered_row = [row[i] for i in cols_to_keep if i < len(row)]
            f.write(OUTPUT_DELIMITER.join(filtered_row) + '\n')

    # Calculate stats for logging
    dropped = max_cols - len(cols_to_keep)
    drop_msg = f"{LogStyle.RED}Dropped {dropped}{LogStyle.RESET}" if dropped > 0 else f"{LogStyle.YELLOW}Dropped 0{LogStyle.RESET}"

    log_success(f"Processed {LogStyle.BOLD}'{filepath}'{LogStyle.RESET} "
                f"-> Kept {len(cols_to_keep)}/{max_cols} cols ({drop_msg}) "
                f"-> Saved as {LogStyle.BOLD}'{out_filepath}'{LogStyle.RESET}")

if __name__ == "__main__":
    print(f"\n{LogStyle.BOLD}--- ZERO-COLUMN FILTER SCRIPT ---{LogStyle.RESET}\n")

    # Filter out files that already have the output suffix to avoid double processing
    files_to_process = [f for f in INPUT_FILES if not f.endswith(OUTPUT_SUFFIX)]

    if not files_to_process:
        log_warn("No new CSV files found matching the criteria.")
    else:
        log_info(f"Found {len(files_to_process)} file(s) to process.\n")
        for file in files_to_process:
            process_csv(file)

    print(f"\n{LogStyle.BOLD}--- DONE ---{LogStyle.RESET}\n")
