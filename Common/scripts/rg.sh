#!/bin/bash
#------------------------------------
# Execure rg.tcl and plot data
#------------------------------------
# To use default frame range:
# => echo -ne "\n\n\n" | ./rg.sh

# Calculate RG
vmd -dispdev text -e rg.tcl

# Plot RG data
python3 plot_data.py -t="Radius of Gyration" -o="rg.pdf" "rg.csv"
