set pdb ../../tclForcesFiles/common/ubiquitin.pdb
set psf ../../tclForcesFiles/common/ubiquitin.psf
set out ubiquitin-grid.pdb

mol load psf $psf pdb $pdb

set all [atomselect top all]
$all set beta 1
$all set occupancy [$all get mass]
$all writepdb $out

exit
