#!/bin/bash

# Calculate RMSD
vmd -dispdev text -e rmsd.tcl

# Plot RMSD data
python3 plot_data.py  -t="RMSD Plot" -o="rmsd.pdf" "rmsd.csv"
