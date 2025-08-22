set pdb ../common/ubiquitin_solvate.pdb
set psf ../common/ubiquitin_solvate.psf
set targetPdb ../common/ubiquitin_solvate_backbone.pdb
set selection "protein and backbone"
set targetMark "1.00"

mol load psf $psf pdb $pdb

set all [atomselect top all]
$all set beta 0
$all set occupancy 0

set target [atomselect top $selection]
set masses [$target get mass]

$target set beta $targetMark
$target set occupancy $masses

$all writepdb $targetPdb

exit
