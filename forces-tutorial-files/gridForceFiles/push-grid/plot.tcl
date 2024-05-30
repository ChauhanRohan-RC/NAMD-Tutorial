set psf		../../tclForcesFiles/common/ubiquitin.psf
set dcd		push-grid.dcd

package require multiplot

set molid [mol load psf $psf dcd $dcd]

set x {}
set y {}
set all [atomselect $molid all]
for { set i 0 } { $i < [molinfo top get numframes] } { incr i } {
    animate goto $i
    set com [measure center $all weight mass]
    lappend x $i
    lappend y [lindex $com 0]
}

set plothandle [multiplot -x $x -y $y -ylabel "Center of mass" -lines -marker point -plot]

$all delete
#mol delete $molid
