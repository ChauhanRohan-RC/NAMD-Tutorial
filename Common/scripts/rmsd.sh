#!/bin/bash
#------------------------------------
# Execure rmsd.tcl and plot data
#------------------------------------
# To use default frame range:
# => echo -ne "\n\n\n" | ./rmsd.sh

# Calculate RMSD
vmd -dispdev text -e rmsd.tcl

# Plot RMSD data
python3 plot_data.py  -t="RMSD Plot" -o="rmsd.pdf" "rmsd.csv"
