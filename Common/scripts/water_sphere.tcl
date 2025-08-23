##############################################################################
### Script to immerse a molecule inside a water sphere						##
##############################################################################

# 1. Copy script in the same directory as input .pdb and .psf files
# 2. INPUT: set molecule name (name of the input .psf and .pdb files) and padding
# 3. run with "vmd -dispdev text -e water_box.tcl"
# 4. OUTPUT: Generates "${molname}_ws" .pdb and .psf files, and .comments.txt file

## KNOWN BUGS:
# -> AUTO-IONIZE: adds ions with residue id's starting from 1.
#    -> DO NOT CHANGE them manually in PDB, will create conflicts with PSF
#    -> To avoid conflicts, use specific selections "nucleic and resid ..." 

# ----------------------------------
# INPUT PARAMS
# ----------------------------------
set molname 			"dna";	# TODO: set molecule name (name of .psf and .pdb files)
set padding 			15;		# TODO: padding (in Å)
set msm_grid_pad		2;		# Multistep Summation Method (MSM) Grid extra padding (in Å) (for safety)

set boundary		2.4;		#  [in Å] min distance to keep between solute and water. Default: 2.4 Å

# ----------------------------------
# OUTPUT PARAMS
# ----------------------------------
set output_name 		"${molname}_ws";		# name of output .psf and .pdb files (defaults to adding "ws" suffix)
set meta_file_name "${molname}_ws.comments.txt";		# name of the meta information file

# Auto-IONIZE (to nuetralize the system) ---------------
set auto_ionize		on;		
set cation			"SOD";		# SOD (Na+), MG (Mg+2), POT (K+), CES (Cs+), CAL (Ca+2), ZN2 (Zn+2)  
set anion			"CLA";		# CLA (Cl-)
set ion_seg_name		"ION";		# New Segment name for ions 

# BETA and OCCUPANCY -------
set set_beta		on;		# TODO: [on/off] whether to set beta values of all atoms
set beta_value	0;		# beta value to set (if $set_beta is ON)

set set_occupancy		on;		# TODO: [on/off] whether to set occupancy values of all atoms
set occupancy_value		0;		# occupancy value to set (if $set_occupancy is ON)




# --------------------------------------------
# SOLVATE
# --------------------------------------------
set mol_id [mol new "${molname}.psf"]
mol addfile "${molname}.pdb" molid $mol_id

### Determine the center of mass of the molecule and store the coordinates
set cen [measure center [atomselect $mol_id all] weight mass]
set x1 [lindex $cen 0]
set y1 [lindex $cen 1]
set z1 [lindex $cen 2]
set max_sq 0

### Determine the distance of the farthest atom from the center of mass
foreach atom [[atomselect $mol_id all] get index] {
	set pos [lindex [[atomselect $mol_id "index $atom"] get {x y z}] 0]
	set x2 [lindex $pos 0]
	set y2 [lindex $pos 1]
	set z2 [lindex $pos 2]
	
	set dist_sq [expr (($x2-$x1)**2) + (($y2-$y1)**2) + (($z2-$z1)**2)]
	if {$dist_sq > $max_sq} { 
		set max_sq $dist_sq 
	}
}

set max [expr sqrt($max_sq)]
set final_r [expr $max + $padding]
set final_r_sq [expr $final_r * $final_r]

mol delete $mol_id

### Solvate the molecule in a water Cube with body_diagonal = sqrt(3) * ($max + $padding) so that it can 
# fully embed the sphere we want
# Hence, effective padding (padding from farthest atom) = body_diagonal - $max = 0.732 * $max + 1.732 * $padding
set pad_eff [expr 0.732 * $max + 1.732 * $padding]

set temp_out "${output_name}.temp"

package require solvate
solvate "${molname}.psf" "${molname}.pdb" -t $pad_eff -b $boundary -o $temp_out

resetpsf
package require psfgen
set mol_id [mol new "${temp_out}.psf"]
mol addfile "${temp_out}.pdb" molid $mol_id
readpsf "${temp_out}.psf"
coordpdb "${temp_out}.pdb"

### Determine which water molecules need to be deleted and use a for loop to delete them
#set wat [atomselect $mol_id "same residue as {water and ((x-$x1)*(x-$x1) + (y-$y1)*(y-$y1) + (z-$z1)*(z-$z1))<($max*$max)}"]
set del [atomselect $mol_id "water and not same residue as { water and ((x-$x1)**2 + (y-$y1)**2 + (z-$z1)**2) < $final_r_sq }"]

set seg [$del get segid]
set res [$del get resid]
set name [$del get name]
for {set i 0} {$i < [llength $seg]} {incr i} {
	delatom [lindex $seg $i] [lindex $res $i] [lindex $name $i] 
}

writepsf "${output_name}.psf"
writepdb "${output_name}.pdb"

# Cleaning up
mol delete $mol_id
resetpsf
rm "${temp_out}.psf" "${temp_out}.pdb" "${temp_out}.log"


# ---------------------------------
# Ionization
# ---------------------------------
if {[string trim $auto_ionize] eq "on"} {
	# Copying solvated .psf and .pdb to temp files (to be used for input here)
	set temp_file_name	"${molname}_ws_solvated_temp"
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
set mol_id [mol new "${output_name}.psf"]
mol addfile "${output_name}.pdb" molid $mol_id

set everyone [atomselect $mol_id all]

# Geometric Center of All Atoms
set cen_geo [measure center $everyone]

# Center of Mass of All Atoms [RECOMMENDED]
set cen_mass [measure center $everyone weight mass]

# Minimum and Maximum Coordinates of All Atoms
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
puts "# ==================================================="
log "# ============== SOLVATE: WATER SPHERE =============="
log "# -> INPUT molname: \"${molname}\" (used for input .psf and .pdb files)"
log "# -> INPUT padding: ${padding} Å"
log "# -> INPUT MSM padding: ${msm_grid_pad} Å"
log "# -> INPUT boundary: ${boundary} Å"
log "# ----------------------------------------------------"
log "# -> INPUT Auto-IONIZE: $auto_ionize | Cation: \"$cation\" | Anion: \"$anion\" | Ion segment name: \"$ion_seg_name\""
log "# -> INPUT Set BETA: $set_beta | BETA Value: $beta_value"
log "# -> INPUT Set OCCUPANCY: $set_occupancy | OCCUPANCY Value: $occupancy_value"
log "# ----------------------------------------------------"
log "# -> OUTPUT solvated: \"${output_name}.psf\", \"${output_name}.pdb\""
log "# NOTE: All dimensions are in Å"
log "# NOTE: Sphere COM should be used in most cases"
log "# ---------------------------------------------------------------"
log "#"
log "# -> Distance of farthest Atom (excluding Water) from COM (Å): $max"
log "# -> SPHERE Padding: ${padding} Å"
log "# -----------------------------"
log "# -> SPHERE Geometric Center: $cen_geo"
log "# -> SPHERE Center of Mass (COM): $cen_mass"
log "# -> SPHERE Radius (Å): $final_r"
log "# -> SPHERE Min Coordinate: $min_vec"
log "# -> SPHERE Max Coordinate: $max_vec"
log "#"
log "# ---------------------------------------------------"
log "# -> COPY-PASTE For NAMD config file "
log "# ---------------------------------------------------"
log "## -> Multilevel Summation Method (MSM) Grid Size with Padding ${msm_grid_pad} Å"
log "MSMxmin [expr floor($xmin) - $msm_grid_pad]"
log "MSMymin [expr floor($ymin) - $msm_grid_pad]"
log "MSMzmin [expr floor($zmin) - $msm_grid_pad]"
log "MSMxmax [expr floor($xmax) + $msm_grid_pad]"
log "MSMymax [expr floor($ymax) + $msm_grid_pad]"
log "MSMzmax [expr floor($zmax) + $msm_grid_pad]"
log "# -------------------------------------------"
log "## -> Spherical Boundary Condition (using COM) (in Å)"
log "sphericalBCcenter   $cen_mass"
log "sphericalBCr1       $final_r"
log "# ---------------------------------------------------------------"
puts ""

# gc
mol delete $mol_id
exit
