wrapmode cell

# Two first agruments of calcforces are automatically forwarded 
# to it by NAMD. The other 4 arguments match the list of 4 values 
# from command tclBCArgs.

proc calcforces {step unique Rstart Rtarget Rrate K} {

  global bubbleCenter ;# defined in tclBCScript{ ... } 

  # increase R, starting from $Rstart, by $Rrate at each step,
  # until it reaches $Rtarget; then keep it constant
  
  set R [expr $Rstart + $Rrate * $step]
  if { $R > $Rtarget } { set R $Rtarget }

  # let only the main processor print the output
  
  if { $unique } {
    print "step $step, bubble radius = $R"
  }

  # get the components of the bubble center vector

  foreach { x0 y0 z0 } $bubbleCenter { break }

  # pick atoms of the given patch one by one

  while {[nextatom]} { 
    set rvec [getcoord] ;# get the atom's coordinates
	
    foreach { x y z } $rvec { break } ;# get components of the vector
	
    # find the distance between the atom and the bubble center
	# (long lines can be broken by a backslash and continued 
	# on the next line)
	
    set rho [expr sqrt(($x-$x0)*($x-$x0) + ($y-$y0)*($y-$y0) + \
      ($z-$z0)*($z-$z0))]
	  
	# if the atom is inside the sphere, push it away radially
	
    if { $rho < $R } {
      set forceX [expr $K * ($x-$x0) / $rho]
      set forceY [expr $K * ($y-$y0) / $rho]
      set forceZ [expr $K * ($z-$z0) / $rho]
      addforce "$forceX $forceY $forceZ"
    }
  }
}


