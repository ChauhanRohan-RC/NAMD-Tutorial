#######################################################################################
### Script to immerse a molecule inside a water box  							 	 ##
#######################################################################################

# 1. Copy script in the same directory as input .pdb and .psf files
# 2. INPUT: set molecule name (name of the input .psf and .pdb files) and padding
# 3. run with "vmd -dispdev text -e water_box.tcl"
# 4. OUTPUT: Generates "${molname}_wb" .pdb and .psf files, and .comments.txt file

## KNOWN BUGS:
# -> AUTO-IONIZE: adds ions with residue id's starting from 1.
#    -> DO NOT CHANGE them manually in PDB, will create conflicts with PSF
#    -> To avoid conflicts, use specific selections "nucleic and resid ..." 

# ----------------------------------
# INPUT PARAMS
# ----------------------------------
set molname		"dna";	# TODO: set molecule name (name of the input .psf and .pdb files)
set padding		15;		# TODO: padding [in Å]. Must be sufficiently high so that system do not interact with its own image in the next cell if PBC is used
set boundary		2.4;		#  [in Å] min distance to keep between solute and water. Default: 2.4 Å

set cubic_box		on;		# TODO: Make a Cubic box of water (with given min padding)
set use_cen_mass		off;		# Use Center of mass (COM) to center solute in water cube, else use Geometric center

# ----------------------------------
# OUTPUT PARAMS
# ----------------------------------
set output_name 		"${molname}_wb";		# name of output .psf and .pdb files (defaults to adding "wb" suffix)
set meta_file_name "${molname}_wb.comments.txt";		# name of the meta information file

# Auto-IONIZE (to nuetralize the system) ---------------
set auto_ionize		on;		
set cation			"SOD";		# SOD (Na+), MG (Mg+2), POT (K+), CES (Cs+), CAL (Ca+2), ZN2 (Zn+2)  
set anion			"CLA";		# CLA (Cl-)
set ion_seg_name		"ION";		# New Segment name for ions 

### BETA and OCCUPANCY -------
set set_beta		on;		# TODO: [on/off] whether to set beta values of all atoms
set beta_value	0;		# beta value to set (if $set_beta is ON)

set set_occupancy		on;		# TODO: [on/off] whether to set occupancy values of all atoms
set occupancy_value		0;		# occupancy value to set (if $set_occupancy is ON)




# ----------------------------------
## PERFORM SOLVATION
# ----------------------------------
package require solvate

# Cubic water box
if {[string trim $cubic_box] eq "on"} {

	# Measure minmax
	set mol_id [mol new "${molname}.psf"]
	mol addfile "${molname}.pdb" molid $mol_id
	
	set everyone [atomselect $mol_id all];
	
	# Min and Max coordinates
	set min_max [measure minmax $everyone];
	set min_vec [lindex $min_max 0];
	set max_vec [lindex $min_max 1];

	set xmin [lindex $min_vec 0]
	set ymin [lindex $min_vec 1]
	set zmin [lindex $min_vec 2]
	set xmax [lindex $max_vec 0]
	set ymax [lindex $max_vec 1]
	set zmax [lindex $max_vec 2]
	
	# length in each dir
	set xlen [expr $xmax - $xmin];
	set ylen [expr $ymax - $ymin];
	set zlen [expr $zmax - $zmin];
	
	# apply padding to max length
	set max_len [expr max($xlen, $ylen, $zlen)];
	set cube_len [expr $max_len + (2 * $padding)];
	set cube_len_half [expr $cube_len / 2];
	
	# Center
	if {[string trim $use_cen_mass] eq "on"} {
		set cen [measure center $everyone weight mass];	# center of mass (COM)
	} else {
		set cen [measure center $everyone];				# Geometric center
	}
	
	# Water Cube Coordinates
	set cube_min [list [expr [lindex $cen 0] - $cube_len_half] [expr [lindex $cen 1] - $cube_len_half] [expr [lindex $cen 2] - $cube_len_half]];
	set cube_max [list [expr [lindex $cen 0] + $cube_len_half] [expr [lindex $cen 1] + $cube_len_half] [expr [lindex $cen 2] + $cube_len_half]];
	
	set cube_minmax [list $cube_min $cube_max];
	
	# gc
	mol delete $mol_id;
	
	# Solvating
	solvate "${molname}.psf" "${molname}.pdb" -minmax $cube_minmax -t 0 -b $boundary -o $output_name;
} else {

	solvate "${molname}.psf" "${molname}.pdb" -t $padding -b $boundary -o $output_name;
}


# ---------------------------------
# Ionization
# ---------------------------------
if {[string trim $auto_ionize] eq "on"} {
	# Copying solvated .psf and .pdb to temp files (to be used for input here)
	set temp_file_name	"${molname}_wb_solvated_temp"
	cp "${output_name}.psf" "${temp_file_name}.psf"
	cp "${output_name}.pdb" "${temp_file_name}.pdb"

	puts "# -> LOG: Performing Ionization to neutralize the system...";
	package require autoionize
	autoionize -psf "${temp_file_name}.psf" -pdb "${temp_file_name}.pdb" -neutralize -cation $cation -anion $anion -seg $ion_seg_name -o $output_name
	
	# Cleaning temp files
	rm "${temp_file_name}.psf" "${temp_file_name}.pdb"
}


# ----------------------------------
# Post-Processing
# ----------------------------------
# Determine the center and minimum-maximum cooridnates of the solvated system (useful for periodic boundary conditions)
set mol_id [mol new "${output_name}.psf"]
mol addfile "${output_name}.pdb" molid $mol_id

set is_pdb_changed	0; # If changed, rewrite pdb file

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

# gc
mol delete $mol_id

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

puts ""
puts ""
puts "# ===================================================="
log "# ================ SOLVATE: WATER BOX ================"
log "# -> INPUT molname: \"${molname}\" (used for input .psf and .pdb files)"
log "# -> INPUT padding: ${padding} Å"
log "# -> INPUT boundary: ${boundary} Å"
log "# -> INPUT cubic_box: ${cubic_box}  |  use_center_mass: ${use_cen_mass}"
log "# ----------------------------------------------------"
log "# -> INPUT Auto-IONIZE: $auto_ionize | Cation: \"$cation\" | Anion: \"$anion\" | Ion segment name: \"$ion_seg_name\""
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

exit
