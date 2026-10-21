#!/bin/bash

# -----------------------------------------------------
# Run with sudo if VMD Lib folder requires it
# -----------------------------------------------------

# TODO: set VMD lib folder
VMD_LIB_DIR="/usr/local/lib/vmd-2.0.0"



# -----------------------------------------------------
if [[ ! -d "$VMD_LIB_DIR" ]]; then
    echo "VMD lib folder not found: ${VMD_LIB_DIR}"
    exit 1
fi

target_dir="${VMD_LIB_DIR}/plugins/noarch/tcl/bigdcd2.0"
echo "Installing bigdcd2.0 to path \"$target_dir\""
echo ""

# copy scripts
mkdir -p "$target_dir"
cp *.tcl "$target_dir"

# file permissions
chmod 644 ${target_dir}/*
