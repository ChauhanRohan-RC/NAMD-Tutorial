# Want to force all atoms
set numatoms 1231

set atoms {}
for { set i 1 } { $i <= $numatoms } { incr i } {
    lappend atoms $i
}

# Tell NAMD the atoms we'll force
foreach atom $atoms {
    addatom $atom
}

# Take acceleration factors from NAMD config file, convert to kcal/(mol*Ang*amu)
set linaccel_namd [vecscale [expr 1.0/418.68] $linaccel]
set angaccel_namd [expr double($angaccel)/418.68]
set PI 3.1415926535898

print "Linear acceleration applied: ($linaccel) Ang*ps^-2"
print "Angular acceleration applied: (0 0 $angaccel) Rad*ps^-2"

proc calcforces { } {
    global atoms numatoms linaccel_namd angaccel_namd
    global PI
    
    # Get coordinates and masses
    loadcoords coords
    loadmasses masses
    
    # First calculate center of mass
    set comsum "0 0 0"
    set totalmass 0
    foreach atom $atoms {
	set comsum [vecadd $comsum [vecscale $masses($atom) $coords($atom)]]
	set totalmass [expr $totalmass + $masses($atom)]
    }
    set com [vecscale [expr 1.0/$totalmass] $comsum]
    print "Center of mass = $com"
    
    foreach atom $atoms {
	# Linear force
	set linforce [vecscale $masses($atom) $linaccel_namd]
	
	# Angular force
	set r [vecsub $coords($atom) $com]
	set x [lindex $r 0]
	set y [lindex $r 1]
	set rho [expr sqrt(pow($x, 2) + pow($y, 2))]
	set phi [expr atan2($y, $x) + $PI/2]
	if { $atom == 1 } {
	    print "atom $atom: phi = $phi"
	}
	set angdir "[expr cos($phi)] [expr sin($phi)] 0.0"
	set angforce [vecscale [expr $masses($atom) * $angaccel_namd * $rho] $angdir]
	
	# Total force
	set force [vecadd $linforce $angforce]
	
	addforce $atom $force
    }
}
