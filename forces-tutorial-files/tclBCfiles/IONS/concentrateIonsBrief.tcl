# we will use the following variable to calculate and print the average 
# number of ions found outside the sphere at each step

set avgNumIons 0

wrapmode cell

##################################################

proc calcforces {step unique} {

  global sphereCenter sphereRadius maxForce forceCoef avgNumIons

  if { $step > 0 && $step % 100 == 0 } { 
	  
    set avgNumIons [expr $avgNumIons / 100.]
    print "Step $step, average number of ions outside the sphere: $avgNumIons"
    set avgNumIons 0
	  
    cleardrops 
	  
  }

  while {[nextatom]} { 

    if { [getmass] < 35 || [getmass] > 40 } {
      dropatom ;# not a K or Cl, forget about it
      continue
    }
	
    # vector between the ion and the sphere's center
	
    set relativePosition [vecsub [getcoord] $sphereCenter] 
    set rho              [veclength $relativePosition]

    if { $rho > $sphereRadius } {
	  addforce [vecscale $relativePosition [expr -$maxForce/$rho]]
	  incr avgNumIons
	} else {
	  dropatom ;# this ion is already inside the sphere
    }
  }
}
