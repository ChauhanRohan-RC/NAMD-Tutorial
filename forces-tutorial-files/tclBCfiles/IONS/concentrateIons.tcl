set ionList {} ;# a list of atom IDs of the ions
set atomID 1   ;# in NAMD, atom IDs are 1-based numbers

set inStream [open $pdbSource r]	;# opening the PDB file

# reading the PDB file line after line:
foreach line [split [read $inStream] \n] {

  # symbols 0-3 should be ATOM or HETA 
  set string1 [string range $line 0 3] 
  
  # symbols 17-20 contain the residue name
  set string2 [string range $line 17 20]  
  
  set string2 [string trim $string2] ;# trimming the residue name
   
  if { [string equal $string1 {ATOM}] || \
	   [string equal $string1 {HETA}] } { 
		   
    # so it's a valid atom entry; let's see if this is an ion
	
    if { [string equal $string2 {CLA}] || \
	     [string equal $string2 {POT}] } {
		 
      # yes it's an ion; append its index to the list
      lappend ionList $atomID
    }
	
    # increase the atomID counter by 1 (this is the default increment)
    incr atomID 
  }
}
close $inStream

#####################################################################

# will use the following variable to calculate and print the average 
# number of ions found outside the sphere at each step

wrapmode cell

# in this case the command tclBCArgs {} located in the NAMD config 
# file doesn't pass any arguments to calcforces, therefore only 
# step and unique will be listed below:

proc calcforces {step unique} {

  # list all global variables set outside this procedure 
  # (except those we won't need, like pdbSource)

  global sphereCenter sphereRadius maxForce ionList 

  # find the components of the sphere's center
  
  foreach { x0 y0 z0 } $sphereCenter { break }
  
  while {[nextatom]} { 
  
  	set atomid [getid] ;# get the ID of the current atom
	
    # check if this ID is listed in ionList; if it's found, 
    # lsearch will return the position of the search pattern 
    # in the list (a number >= 0), otherwise -1 is returned
	
	if { [lsearch $ionList $atomid] >= 0 } {
  
	  set rvec [getcoord] ;# get the vector with ion's coordinates
      foreach { x y z } $rvec { break } ;# get the vector's components
	
      # find the distance between the ion and the sphere's center

      set rho [expr sqrt(($x-$x0)*($x-$x0) + ($y-$y0)*($y-$y0) + \
        ($z-$z0)*($z-$z0))]
		
      # Apply same force $maxForce to each ion if it's outside the
      # sphere. The components of the force vector are chosen so 
      # that the vector is directed toward the sphere's center, 
      # and the vector's norm is $maxForce.
	  
      if { $rho > $sphereRadius } {
	  
        set forceX [expr -$maxForce * ($x-$x0) / $rho]
        set forceY [expr -$maxForce * ($y-$y0) / $rho]
        set forceZ [expr -$maxForce * ($z-$z0) / $rho]
		
      } else {
	  
      # If the ion in already inside the sphere, scale the force 
      # by a factor of $rho/$sphereRadius, so that the force 
      # decreases from $maxForce to 0 as the ion approaches the 
      # sphere's center:
	  
        set forceX [expr -$maxForce * ($x-$x0) / $sphereRadius]
        set forceY [expr -$maxForce * ($y-$y0) / $sphereRadius]
        set forceZ [expr -$maxForce * ($z-$z0) / $sphereRadius]
		
      }

      # Finally, applying the force calculated above 
      # to the current atom

      addforce "$forceX $forceY $forceZ"
    }
  }
}
