#!/bin/bash

# Calculate RG
vmd -dispdev text -e rg.tcl

# Plot RG data
python3 plot_data.py -t="Radius of Gyration" -o="rg.pdf" "rg.csv"
