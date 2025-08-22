#######################################################################################
### Script to immerse a molecule inside a water box  							 	 ##
#######################################################################################

# 1. Copy script in the same directory as input .pdb and .psf files
# 2. INPUT: set molecule name (name of the input .psf and .pdb files) and padding
# 3. run with "vmd -dispdev text -e water_box.tcl"
# 4. OUTPUT: Generates "${molname}_wb" .pdb and .psf files, and .comments.txt file

# INPUT PARAMS ---------------------------
set molname "dna";	# TODO: set molecule name (name of the input .psf and .pdb files)
set padding 5;		# TODO: padding in Angrostom


# OUTPUT PARAMS --------------------------
set output_name "${molname}_wb";		# name of output .psf and .pdb files (defaults to adding "wb" suffix)

set meta_file_name "${molname}_wb.comments.txt";		# name of the meta information file

### BETA and OCCUPANCY -------
set set_beta		on;		# TODO: [on/off] whether to set beta values of all atoms
set beta_value	0;		# beta value to set (if $set_beta is ON)

set set_occupancy		on;		# TODO: [on/off] whether to set occupancy values of all atoms
set occupancy_value		0;		# occupancy value to set (if $set_occupancy is ON)




# ----------------------------------------------------
## PERFORM SOLVATION
# ----------------------------
package require solvate
solvate "${molname}.psf" "${molname}.pdb" -t $padding -o $output_name


# -------------------------------
# META INFO
# -------------------------------

# open meta file for writing
set meta_file [open $meta_file_name w]

proc log2meta { msg } { 
	global meta_file;
	puts $meta_file $msg;	# to meta file
}

proc log { msg } { 
	global meta_file;
	puts $msg; flush stdout;		# to stdout
	log2meta $msg;					# to meat file
}

# ------------ FINAL LOAD ------------------------
# Determine the center and minimum-maximum cooridnates of the solvated system (useful for periodic boundary conditions)
set mol_id [mol new "${output_name}.psf"]
mol addfile "${output_name}.pdb" molid $mol_id

set is_pdb_changed 0; # If changed, rewrite pdb file

set everyone [atomselect $mol_id all]

# Geometric Center of all atoms (system + water)
set cen_geo [measure center $everyone]

# Center of Mass of all atoms (system + water)
set cen_mass [measure center $everyone weight mass]

# Minimum and Maximum Coordinates of water box
set min_max [measure minmax $everyone]
set min_vec [lindex $min_max 0]
set max_vec [lindex $min_max 1]

set xmin [lindex $min_vec 0]
set ymin [lindex $min_vec 1]
set zmin [lindex $min_vec 2]
set xmax [lindex $max_vec 0]
set ymax [lindex $max_vec 1]
set zmax [lindex $max_vec 2]

# Setting BETA and OCCUPANCY --------------
puts ""
if {[string trim $set_beta] eq "on"} {
	puts "# -> LOG: SETTING BETA = $beta_value";
	$everyone set beta $beta_value;
	set is_pdb_changed 1;
}

if {[string trim $set_occupancy] eq "on"} {
	puts "# -> LOG: SETTING OCCCUPANCY = $occupancy_value";
	$everyone set occupancy $occupancy_value;
	set is_pdb_changed 1;
}


# IF PDB is changed, write it again
if { $is_pdb_changed == 1 } {
	puts "# -> LOG: Rewriting PDB file ${output_name}.pdb";
	puts ""
	$everyone writepdb "${output_name}.pdb";
}

puts ""
puts ""
puts "# ===================================================="
log "# ================ SOLVATE: WATER BOX ================"
log "# -> INPUT molname: \"${molname}\" (used for input .psf and .pdb files)"
log "# -> INPUT padding: ${padding} Å"
log "# ----------------------------------------------------"
log "# -> INPUT Set BETA: $set_beta | BETA Value: $beta_value"
log "# -> INPUT Set OCCUPANCY: $set_occupancy | OCCUPANCY Value: $occupancy_value"
log "# ----------------------------------------------------"
log "# -> OUTPUT solvated: \"${output_name}.psf\", \"${output_name}.pdb\""
log "# NOTE: All dimensions are in Å"
log "#-----------------------------------------------------"
log "#"
log "# -> Geometric Center of All Atoms: $cen_geo"
log "# -> Center of Mass (COM) of All Atoms: $cen_mass"
log "# -> Box Min Coordinate: $min_vec"
log "# -> Box Max Coordinate: $max_vec"
log "#"
log "# ---------------------------------------------------"
log "# -> COPY-PASTE For NAMD config file (Periodic BC): "
log "# ---------------------------------------------------"
log "cellBasisVector1 [ expr $xmax - $xmin ] 0 0 "
log "cellBasisVector2 0 [ expr $ymax - $ymin ] 0 "
log "cellBasisVector3 0 0 [ expr $zmax - $zmin ] "
log "cellOrigin [ expr ($xmax + $xmin) / 2 ] [ expr ($ymax + $ymin) / 2 ] [ expr ($zmax + $zmin) / 2 ] "
log "# ---------------------------------------------------"
puts ""

# gc
mol delete $mol_id
exit
