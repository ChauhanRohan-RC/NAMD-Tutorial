set atomList {} ;# a list of atom IDs of water oxigens
set atomID 1   ;# in NAMD, atom IDs are 1-based numbers

set inStream [open $pdbSource r]	;# opening the PDB file

# reading the PDB file line after line:
foreach line [split [read $inStream] \n] {

  set string1 [string range $line 0 3]    ;# symbols 0-3 should be ATOM or HETA 
  set string2 [string range $line 13 16]  ;# symbols 13-16 contain the atom's name
  set string2 [string trim $string2]      ;# trimming the atom's name
   
  if { [string equal $string1 {ATOM}] || \
	   [string equal $string1 {HETA}] } { 
		   
    # so it's a valid atom entry; let's see if this is a water's oxigen
	
    if { [string equal $string2 {OH2}] } {

      lappend atomList $atomID ;# yes it is; append its index to the list
    }
	incr atomID ;# increase the atomID counter by 1 (this is the default increment)
  }
}
close $inStream 

puts "[llength $atomList] water oxigens atoms found"

set totalTorque 0.

wrapmode cell

###########################################################################

proc calcforces {step unique} {

  global atomList shearStress zLo zHi dz TOL force

  # the values of global variables is stored between calls to a procedure
  global totalTorque 

  if { $step % 100 == 0 } { 
  
    if { $totalTorque != 0. } { 
      print "Step $step, total torque applied: $totalTorque" 
	  set totalTorque 0.
	}
    cleardrops 
  }

  while {[nextatom]} { 

    # check if this atom should be considered at all

	if { [lsearch $atomList [getid]] == -1 } {
	  dropatom
	  continue
    }

    # now check if it's within bounds zLo to zHi

    foreach { x y z } [getcoord] { break } ;# get coordinates
	
        if { $z >= [expr {$zLo-$dz}] && $z <= [expr {$zLo+$dz}] } {
	
	  addforce "[expr {-$force}] 0.0 0.0"
	  set totalTorque [expr {$totalTorque - $force * $z}]

	} elseif { $z >= [expr {$zHi-$dz}] && $z <= [expr {$zHi+$dz}] } {
	
	  addforce "$force 0.0 0.0"
	  set totalTorque [expr {$totalTorque + $force * $z}]

	} elseif { ( $z >= [expr {$zLo-$TOL}] && $z <= [expr {$zLo+$TOL}] ) || \
               ( $z >= [expr {$zHi-$TOL}] && $z <= [expr {$zHi+$TOL}] ) } {

      continue ;# keep an eye on them, they may come closer next time

    } else {
	
      dropatom ;# this atom is too far, forget about it for a while
	  
    }
  }
}
