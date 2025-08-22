#!/bin/bash

# Calculate RMSF
vmd -dispdev text -e rmsf.tcl

# Plot RMSF data
python3 plot_data.py -t="RMSF Plot" -o="rmsf.pdf" "rmsf.csv"
