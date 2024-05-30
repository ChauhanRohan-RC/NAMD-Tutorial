set pdb graphene-ion.pdb
set psf graphene-ion.psf
set out graphene-grid.pdb

mol load psf $psf pdb $pdb

set all [atomselect top all]
$all set beta 0
$all set occupancy 0

set sel [atomselect top "noh"]
$sel set beta 1
$sel set occupancy 1
$all writepdb $out

exit
