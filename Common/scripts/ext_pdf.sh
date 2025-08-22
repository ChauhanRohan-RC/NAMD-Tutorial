#!/bin/bash
#
# Script to calculate and plot Extension PDF
# TODO
# Set INPUT and OUTPUT params in "distance.tcl" and "ext_pdf.py"
#

# 1. Calculate Extension
vmd -dispdev text -e distance.tcl

# 2. Plot Ext vs Frame (time) and Extension distribution
python3 ext_pdf.py
