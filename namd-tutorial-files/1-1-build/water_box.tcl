### Script to generate a solvation box around a molecule

# TODO: molecule name (also the input file name)
set molname ubq

# TODO: padding in Angrostom
set padding 5

package require solvate
solvate ${molname}.psf ${molname}.pdb -t $padding -o ${molname}_wb

# Determine the center and minimum-maximum cooridnates of the solvated system (useful for periodic boundary conditions)
mol new ${molname}_wb.psf
mol addfile ${molname}_wb.pdb

set everyone [atomselect top all]

# Geometric Center of water box
set cen [measure center $everyone]

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

puts "---------------------------------------------------"
puts "Water Box Geometric Center: $cen"
puts "Water Box Min: $min_vec"
puts "Water Box Max: $max_vec"
puts "---------------------------------------------------"
puts "Copy-Paste For NAMD config file: "
puts "cellBasisVector1 [ expr $xmax - $xmin ] 0 0 "
puts "cellBasisVector2 0 [ expr $ymax - $ymin ] 0 "
puts "cellBasisVector3 0 0 [ expr $zmax - $zmin ] "
puts "cellOrigin [ expr ($xmax + $xmin) / 2 ] [ expr ($ymax + $ymin) / 2 ] [ expr ($zmax + $zmin) / 2 ] "
puts "-------------------------------------------------------"
