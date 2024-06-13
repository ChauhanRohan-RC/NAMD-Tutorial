##############################################################################
### Script to immerse a molecule inside a water sphere						##
##############################################################################

# 1. Set molecule name (name of the input .psf and .pdb files)
# 2. Make sure the script is in the same directory as input .pdb and .psf files
# 2. run with "vmd -dispdev text -e water_sphere.tcl"

set molname ubq;	# TODO: set molecule name (name of .psf and .pdb files)
set padding 10;		# TODO: padding in Angstrom

set output_name ${molname}_ws;	# name of output .psf and .pdb files (defaults to adding "ws" suffix)
set msm_grid_pad 2;				# MSM Grid extra padding in Å (for safety)

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
solvate "${molname}.psf" "${molname}.pdb" -t $pad_eff -o $temp_out

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

mol delete $mol_id


### =================== FINAL LOAD  ======================
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

puts ""
puts ""
puts "===================================="
puts "Output Files: ${output_name}.psf   ${output_name}.pdb"
puts "NOTE: All dimensions are in Å"
puts "NOTE: Sphere COM should be used in most cases"
puts "---------------------------------------------------------------"
puts "Distance of farthest Atom (excluding Water) from COM (Å): $max"
puts "Sphere Padding (Å): $padding"
puts "---------------------------------------------------------------"
puts "Sphere Geometric Center: $cen_geo"
puts "Sphere Center of Mass (COM): $cen_mass"
puts "Sphere Radius (Å): $final_r"
puts "Sphere Min Coordinate: $min_vec"
puts "Sphere Max Coordinate: $max_vec"
puts "---------------------------------------------------------------"
puts "============ Copy-Paste For NAMD config file ================"
puts "## Multilevel Summation Method (MSM) Grid Size with Padding ${msm_grid_pad} Å"
puts "MSMxmin [expr floor($xmin) - $msm_grid_pad]"
puts "MSMymin [expr floor($ymin) - $msm_grid_pad]"
puts "MSMzmin [expr floor($zmin) - $msm_grid_pad]"
puts "MSMxmax [expr floor($xmax) + $msm_grid_pad]"
puts "MSMymax [expr floor($ymax) + $msm_grid_pad]"
puts "MSMzmax [expr floor($zmax) + $msm_grid_pad]"
puts ""
puts "## Spherical Boundary Condition (using COM) (in Å)"
puts "sphericalBCcenter   $cen_mass"
puts "sphericalBCr1       $final_r"
puts "---------------------------------------------------------------"

# gc
mol delete $mol_id
exit
