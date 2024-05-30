#######################################################################################
### Script to immerse a molecule inside a water box  							 	 ##
#######################################################################################

# 1. Set molecule name (name of the input .psf and .pdb files)
# 2. Make sure the script is in the same directory as input .pdb and .psf files
# 2. run with "vmd -dispdev text -e water_box.tcl"

set molname ubq;	# TODO: set molecule name (name of the input .psf and .pdb files)
set padding 5;		# TODO: padding in Angrostom

set output_name ${molname}_wb;		# name of output .psf and .pdb files (defaults to adding "wb" suffix)

package require solvate
solvate ${molname}.psf ${molname}.pdb -t $padding -o $output_name

# Determine the center and minimum-maximum cooridnates of the solvated system (useful for periodic boundary conditions)
mol new ${output_name}.psf
mol addfile ${output_name}.pdb

set everyone [atomselect top all]

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

puts ""
puts ""
puts "============================"
puts "Output Files: ${output_name}.psf   ${output_name}.pdb"
puts "NOTE: All dimensions are in Å"
puts "---------------------------------------------------"
puts "Geometric Center of All Atoms: $cen_geo"
puts "Center of Mass (COM) of All Atoms: $cen_mass"
puts "Box Min Coordinate: $min_vec"
puts "Box Max Coordinate: $max_vec"
puts "---------------------------------------------------"
puts "Copy-Paste For NAMD config file: "
puts "cellBasisVector1 [ expr $xmax - $xmin ] 0 0 "
puts "cellBasisVector2 0 [ expr $ymax - $ymin ] 0 "
puts "cellBasisVector3 0 0 [ expr $zmax - $zmin ] "
puts "cellOrigin [ expr ($xmax + $xmin) / 2 ] [ expr ($ymax + $ymin) / 2 ] [ expr ($zmax + $zmin) / 2 ] "
puts "-------------------------------------------------------"
