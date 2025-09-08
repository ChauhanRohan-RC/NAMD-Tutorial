#!/bin/bash
#------------------------------------
# Execure rmsf.tcl and plot data
#------------------------------------
# To use default frame range:
# => echo -ne "\n\n\n" | ./rmsf.sh

# Calculate RMSF
vmd -dispdev text -e rmsf.tcl

# Plot RMSF data
python3 plot_data.py -t="RMSF Plot" -o="rmsf.pdf" "rmsf.csv"
