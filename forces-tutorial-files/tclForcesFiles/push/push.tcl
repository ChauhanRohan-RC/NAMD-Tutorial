# Want to force all atoms
set numatoms 1231

set atoms {}
for { set i 1 } { $i <= $numatoms } { incr i } {
    lappend atoms $i
}

foreach atom $atoms {
    addatom $atom
}

# Convert input to NAMD units: kcal/(mol*Ang*amu)
set linaccel_namd [vecscale [expr 1.0/418.68] $linaccel]

print "Linear acceleration applied: ($linaccel) Ang*ps^-2"

proc calcforces { } {
    global atoms numatoms linaccel_namd
    
    loadcoords coords
    loadmasses masses
    
    set comsum "0.0 0.0 0.0"
    set totalmass 0.0
    foreach atom $atoms {
        # Take force vector from NAMD config file
        set force [vecscale $masses($atom) $linaccel_namd]
        addforce $atom $force
        set tmp [vecscale $masses($atom) $coords($atom)]
        set comsum [vecadd $comsum $tmp]
        set totalmass [expr $totalmass + $masses($atom)]
    }
    set invmass [expr 1.0/$totalmass]
    print "Center of mass = [vecscale $invmass $comsum]"
}
