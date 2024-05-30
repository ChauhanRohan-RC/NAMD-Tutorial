set targetMark    "1.00"

set targets {}
set masses {}

set inStream [open $targetAtomPdb r]
foreach line [split [read $inStream] \n] {
    set string1 [string range $line 0 3]
    set string2 [string range $line 6 10]
    set string3 [string range $line 17 20]	
    set string4 [string range $line 12 15]
    set string5 [string range $line 46 53]
    set string6 [string range $line 72 75]
    set string673 [string range $line 73 75]
    set string7 [string range $line 22 25]
    set string8 [string range $line 62 65]
    set string9 [string range $line 55 60]
    
    
    if { ([string equal $string1 {ATOM}] || \
 	      [string equal $string1 {HETA}] ) && \
 	     [string equal $targetMark $string8] } {	
 	lappend targets "[string trim $string6]\
 			    [string trim $string7] [string trim $string4]"
	lappend masses "[string trim $string9]"
    }
}
close $inStream


# make list of atoms
set atoms {}
foreach target $targets {
    foreach {segname resid atom} $target { break }
    set atomindex [atomid $segname $resid $atom]
    lappend atoms $atomindex
    addatom $atomindex
}

set numatoms [llength $atoms]

if { $numatoms > 0 } {
    set applyforce 1
} else {
    print "WARNING: no target atoms have been detected"
    set applyforce 0
}

# Take force factor from NAMD config file
set linaccel_namd [vecscale [expr 1.0/418.68] $linaccel]
set angaccel_namd [expr double($angaccel)/418.68]
set PI 3.1415926535898

print "Linear acceleration applied: ($linaccel) Ang*ps^-2"
print "Angular acceleration applied: (0 0 $angaccel) Rad*ps^-2"

# Zero out initial forces
set forces {}
foreach index $atoms {
    lappend forces "0.0 0.0 0.0"
}

set forcecount $forcesRecalcFreq
set printcount 0

proc calcforces { } {
    global atoms numatoms masses linaccel_namd angaccel_namd forces
    global applyforce forcesRecalcFreq
    global forcecount printcount
    global PI
    
    if { $applyforce } {
	# Apply forces
	foreach atom $atoms force $forces {
	    addforce $atom $force
	}
	
	# Check to see if we're about to recalculate forces
	# next time step. If so, clearconfig and re-add atoms
	if { $forcecount == [expr $forcesRecalcFreq - 1] } {
	    print "Adding atoms prior to reconfiguring forces at $printcount"
	    foreach atom $atoms {
		addatom $atom
	    }
	}
	
	if { $forcecount == $forcesRecalcFreq } {
	    print "Recalculating forces at $printcount"
	    
	    # Get coordinates
	    loadcoords coords
	    
	    # First calculate center of mass
	    set comsum "0 0 0"
	    set totalmass 0
	    foreach atom $atoms mass $masses {
		set comsum [vecadd $comsum [vecscale $mass $coords($atom)]]
		set totalmass [expr $totalmass + $mass]
	    }
	    set com [vecscale [expr 1.0/$totalmass] $comsum]
	    print "Center of mass = $com"
	    
	    # Now calculate forces
	    set forces {}
	    foreach atom $atoms mass $masses {
		# Linear force
		set linforce [vecscale $mass $linaccel_namd]
		
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
		set angforce [vecscale [expr $mass * $angaccel_namd * $rho] $angdir]
		
		# Total force
		set force [vecadd $linforce $angforce]
		
		lappend forces $force
	    }
	    
	    print "Step ${printcount}: Recalculated [llength $forces] forces"
	    
	    set forcecount 0
	    
	    # Vital to clearconfig after calculating forces, else
	    # coordinates will be retrieved every timestep, erasing
	    # much of the potential speed gains
	    clearconfig
	}
	incr forcecount
    }
    incr printcount
    return
}
