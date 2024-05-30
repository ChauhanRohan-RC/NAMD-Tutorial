set pdb ../../tclBCexamples/IONS/waterbox40-0.2M.pdb
set psf ../../tclBCexamples/IONS/waterbox40-0.2M.psf
set out1 potassium.pdb
set out2 chloride.pdb

mol load psf $psf pdb $pdb

set all [atomselect top all]
$all set beta 0
$all set occupancy 0

set sel [atomselect top "name POT"]
$sel set beta 1
$sel set occupancy 1
$all writepdb $out1

$all set beta 0
$all set occupancy 0
set sel [atomselect top "name CLA"]
$sel set beta 1
$sel set occupancy 1
$all writepdb $out2

exit
