#!/usr/bin/env python3

import math
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.transformations import wrap, center_in_box
from MDAnalysis.transformations.boxdimensions import set_dimensions
from tqdm import tqdm

"""
Script to Wrap MD trajectories with Periodic Boundary Conditions (PBC Wrap)
Equivalent to "pbc wrap" of VMD PBCTools, but very efficient. 

=> Process trajectory frame by frame
=> Supports multiprocessing

Usage: scroll to bottom and Search for todo 
"""


# ==============================================
# HELPER FUNCTIONS
# ==============================================

def find_files(dir_path, prefix, suffix, min_num=None, max_num=None, sort_natural=True, return_abs_path=False):
    """
    Finds files in a directory matching a specific prefix, numerical sequence, and suffix.

    ex. find_files(".", "amyld_wb_eq", ".coor")
    """
    dir_path_obj = Path(dir_path)

    # Uncomment the following lines to enable the directory check
    # (matching your commented Tcl code)
    if not dir_path_obj.is_dir():
        raise FileNotFoundError(f"Directory '{dir_path}' not found.")

    # Escape prefix and suffix so special characters (like '.') are treated literally
    pattern = re.compile(rf"^{re.escape(prefix)}(\d+){re.escape(suffix)}$")
    result_list = []

    # Iterate through all items in the directory
    for f in dir_path_obj.iterdir():
        if f.is_file():
            match = pattern.search(f.name)
            if match:
                # Extract the number and convert to integer (handles removing leading zeros)
                num = int(match.group(1))

                # Check boundaries (if min_num/max_num are provided)
                if (min_num is None or num >= min_num) and \
                        (max_num is None or num <= max_num):

                    if return_abs_path:
                        # f.resolve() gets the absolute, normalized path
                        fpath = str(f.resolve())
                    else:
                        # str(f) returns the path exactly as it was constructed (dir_path/filename)
                        fpath = str(f.name)

                    result_list.append(fpath)

    # Apply natural/dictionary sorting
    if sort_natural:
        def natural_sort_key(s):
            # Splits the string into chunks of text and numbers, converting numbers to ints
            # This mimics Tcl's `lsort -dictionary`
            return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]

        result_list.sort(key=natural_sort_key)

    return result_list


def parse_xsc(xsc_path):
    """
    Parses a NAMD .xsc file and returns MDAnalysis compatible box dimensions
    for any box type (orthogonal, triclinic, truncated octahedron, etc.).
    Returns: [length_A, length_B, length_C, alpha, beta, gamma]
    """
    try:
        with open(xsc_path, "r") as f:
            # Ignore comments and empty lines
            lines = [line for line in f.readlines() if not line.strip().startswith("#")]

            # The last line holds the step number, the 3 vectors, and the origin
            data = lines[-1].split()

            # Extract the 3 basis vectors (skip data[0] which is the step number)
            Ax, Ay, Az = float(data[1]), float(data[2]), float(data[3])
            Bx, By, Bz = float(data[4]), float(data[5]), float(data[6])
            Cx, Cy, Cz = float(data[7]), float(data[8]), float(data[9])

            # 1. Calculate vector magnitudes (Lengths A, B, C)
            len_A = math.sqrt(Ax ** 2 + Ay ** 2 + Az ** 2)
            len_B = math.sqrt(Bx ** 2 + By ** 2 + Bz ** 2)
            len_C = math.sqrt(Cx ** 2 + Cy ** 2 + Cz ** 2)

            # Prevent division by zero if the box is completely collapsed
            if len_A == 0 or len_B == 0 or len_C == 0:
                raise ValueError("One or more box vectors have a length of exactly zero.")

            # 2. Calculate angles in degrees
            # Note: We use max(-1.0, min(1.0, value)) to clamp the dot product.
            # This prevents floating point rounding errors (like 1.0000000000000002)
            # from crashing the math.acos() function.

            # Alpha: angle between B and C
            dot_BC = (Bx * Cx) + (By * Cy) + (Bz * Cz)
            cos_alpha = dot_BC / (len_B * len_C)
            alpha = math.degrees(math.acos(max(-1.0, min(1.0, cos_alpha))))

            # Beta: angle between A and C
            dot_AC = (Ax * Cx) + (Ay * Cy) + (Az * Cz)
            cos_beta = dot_AC / (len_A * len_C)
            beta = math.degrees(math.acos(max(-1.0, min(1.0, cos_beta))))

            # Gamma: angle between A and B
            dot_AB = (Ax * Bx) + (Ay * By) + (Az * Bz)
            cos_gamma = dot_AB / (len_A * len_B)
            gamma = math.degrees(math.acos(max(-1.0, min(1.0, cos_gamma))))

            # MDAnalysis expects a list of exactly 6 floats
            return [len_A, len_B, len_C, alpha, beta, gamma]

    except Exception as e:
        print(f"CRITICAL ERROR: Failed to parse XSC file '{xsc_path}'.")
        print(f"Details: {e}")
        sys.exit(1)


def process_trajectory(file_path_str, box_dims):
    """
    The worker function executed by each CPU core.
    We pass strings and lists (picklable) and build the Universe inside the worker.
    """
    input_path = Path(file_path_str)

    # Generate the output filename while preserving the original extension
    output_name = f"{out_file_prefix}{input_path.stem}{out_file_suffix}{input_path.suffix}"
    _out_dir = input_path.parent
    if out_dir and os.path.isdir(out_dir):
        _out_dir = os.path.abspath(out_dir)
    output_path = os.path.join(_out_dir, output_name)

    try:
        # 1. Initialize Universe for this specific file
        u = mda.Universe(input_psf, str(input_path))
        anchor_group = u.select_atoms(center_atom_selection)

        # 2. Build the transformation workflow dynamically
        workflow = []

        # If it's a non-DCD file, inject the box dimensions first
        if input_path.suffix.lower() != '.dcd':
            if box_dims is None:
                raise ValueError(f"Box dimensions missing for non-DCD file: {input_path.name}")
            workflow.append(set_dimensions(box_dims))

        # Add the standard centering and wrapping
        workflow.extend([
            center_in_box(anchor_group, center=center_method),
            wrap(u.atoms, compound=wrap_compound)
        ])

        u.trajectory.add_transformations(*workflow)

        # 3. Execute and write the new trajectory
        # mda.Writer automatically infers the correct format from output_path.suffix
        traj_slice = u.trajectory
        with mda.Writer(str(output_path), u.atoms.n_atoms) as W:
            # for ts in traj_slice:
            #     W.write(u.atoms)
            for ts in tqdm(traj_slice, total=len(traj_slice), desc=f"[+] {input_path.name}", unit="frames"):
                W.write(u.atoms)

        return f"[+] SUCCESS: {input_path} -> {output_path}"

    except Exception as e:
        return f"[-] FAILED: {input_path} - Error: {str(e)}"


# ==============================================================================
# 3. MAIN EXECUTION BLOCK
# ==============================================================================

def main():
    print("--- Starting MDAnalysis pBC Wrapping Pipeline ---")

    # 1. Precondition Checks
    if not os.path.exists(input_psf):
        print(f"CRITICAL ERROR: PSF file '{input_psf}' not found.")
        sys.exit(1)

    has_non_dcd = any(Path(f).suffix.lower() != '.dcd' for f in input_traj_files)
    box_dimensions = None

    if has_non_dcd:
        print("Notice: Non-DCD files detected in input array.")
        if not input_xsc_file or not os.path.exists(input_xsc_file):
            print(f"CRITICAL ERROR: An XSC file is required for non-DCD files, but '{input_xsc_file}' was not found.")
            sys.exit(1)

        # 2. Parse XSC
        box_dimensions = parse_xsc(input_xsc_file)
        print(f"Parsed Box Dimensions from XSC: {box_dimensions}")
    else:
        print("Notice: Only DCD files detected. Box dimensions will be read directly from trajectories.")

    # 3. Filter valid files
    valid_files = [f for f in input_traj_files if os.path.exists(f)]
    missing_files = set(input_traj_files) - set(valid_files)

    if missing_files:
        print(f"Warning: The following files were not found and will be skipped: {missing_files}")

    if not valid_files:
        print("No valid trajectory files to process. Exiting.")
        sys.exit(0)

    # 4. Multiprocessing Execution
    print(
        f"\nSpinning up {max_parallel_workers if max_parallel_workers else 'all available'} CPU workers for {len(valid_files)} files...")

    worker_count = min(max_parallel_workers, len(valid_files)) if max_parallel_workers else None
    with ProcessPoolExecutor(max_workers=worker_count) as executor:
        results = executor.map(process_trajectory, valid_files, [box_dimensions] * len(valid_files))

        for result in results:
            print(result)

    print("\n--- PBC Wrap Complete! ---")



# ==============================================================================
# INPUT CONFIGURATION
# ==============================================================================

input_psf = "../common/amyld_wb.psf"  # todo: input psf

# todo: input trajectory file(s)
input_traj_files = find_files(".", "amyld_wb_eq", ".dcd", return_abs_path=True)

# todo: [OPTIONAL] Periodic cell dimensions to use with non-dcd files
input_xsc_file = None

out_dir = None        # [OPTIONAL] output directory. Defaults to the folder of each specific input file
out_file_prefix = ""
out_file_suffix = ".wrapped"  # todo: output file suffix

center_atom_selection = "protein and name CA"  # todo: atom selection which will be centered in the box
center_method = "geometry"  # 'geometry' (unweighted) or 'mass' (weighted)
wrap_compound = "residues"  # [atoms, group, residues, segments, fragments] Atoms in this group are transformed together

max_parallel_workers = None  # or None for all cpu cores. Each traj file only uses 1 core


if __name__ == '__main__':
    main()
