#!/bin/tclsh

###################################################################
## Script to calculate Min and Max coordinates in a PDB file     ##
###################################################################

## USAGE
# 1. copy this script to the working dir
# 2. INPUT: PDB file and ATOM_SELECTION
# 3. OUTPUT: (optional) output file
# 2. run with "vmd -dispdev text -e minmax.tcl"

# INPUT ----------------------------------------
set pdb_file			"ubq.pdb";

# Atoms to select. Use "all" for every atom
set atom_selection		"all";	

# OUTPUT (optional) ----------------------------

# Output file. Comment out for no output
set out_file_name 		"${pdb_file}.minmax.out";

#----------------------------------------------

# Adding molecule
set mol_id [mol new $pdb_file]

# ATom selection
set sel [atomselect $mol_id $atom_selection]

# Min and Max coordinates
set min_max [measure minmax $sel]

set min_vec [lindex $min_max 0]
set max_vec [lindex $min_max 1]

set xmin [lindex $min_vec 0]
set ymin [lindex $min_vec 1]
set zmin [lindex $min_vec 2]

set xmax [lindex $max_vec 0]
set ymax [lindex $max_vec 1]
set zmax [lindex $max_vec 2]

# Cell Length
set len_x [expr $xmax - $xmin]
set len_y [expr $ymax - $ymin]
set len_z [expr $zmax - $zmin]

# Cell center from Origin
set cen_x [expr ($xmax + $xmin) / 2]
set cen_y [expr ($ymax + $ymin) / 2]
set cen_z [expr ($zmax + $zmin) / 2]

# If there is "output_file" defined
if {[info exists out_file_name] && [string trim $out_file_name] != "" } {
	set has_out_file 1
	
	# open output file ofor writing
	set out_file [open $out_file_name w]
} else {
	set has_out_file 0
	set out_file ""
}

proc log2file { msg } { 
	global has_out_file;
	global out_file;
	
	if { $has_out_file == 1 } {
		puts $out_file $msg;	# to output file
	}
}

proc log { msg } {
	puts $msg; flush stdout;		# to stdout
	log2file $msg;					# to output file
}

puts ""
log "================= Min-Max Coordinates =================="
log " -> Input PDB: \"${pdb_file}\""
log " -> ATOM Selection: \"${atom_selection}\"  |  Selected ATOM Count: [$sel num]"
log "-----------------------------------------------"
log " Min and Max coordinates of selected atoms: "
log " -> Min Coordinate: ${min_vec}"
log " -> Max Coordinate: ${max_vec}"
log "-----------------------------------------------"
log " Length in each direction: "
log " -> X: ${len_x}"
log " -> Y: ${len_y}"
log " -> Z: ${len_z}"
log "-----------------------------------------------"
log " -> Cell CENTER (From origin): { $cen_x $cen_y $cen_z }"
log "-------------------------------------------------------"
log " Copy/paste for PERIODIC BC in NAMD config File: "
log "cellBasisVector1		${len_x} 0 0; "
log "cellBasisVector2		0 ${len_y} 0; "
log "cellBasisVector3		0 0 ${len_z}; "
log "cellOrigin			$cen_x $cen_y $cen_z; "
puts "====================================================="

if { $has_out_file == 1 } {
	close $out_file;
	puts " FINISHED  |  Output File: ${out_file_name}";
	puts "================================================="
} 

puts ""
exit
