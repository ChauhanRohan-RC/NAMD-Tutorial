wrapmode cell

# this line is moved above calcforces to avoid doing this 
# transformation more than once

foreach { x0 y0 z0 } $bubbleCenter { break }

set TOL 3. ;# distance tolerance parameter

proc calcforces {step unique Rstart Rtarget Rrate K} {

  global x0 y0 z0 TOL

  set R [expr {$Rstart + $Rrate * $step}]
  if { $R > $Rtarget } { set R $Rtarget }

  if { $unique } {
    print "step $step, bubble radius = $R"
  }

  # Atoms found at a distance larger than $TOL from the surface
  # of the bubble, or $RTOL from the center of the bubble, will
  # be ignored (dropped) for the rest of a 100-step cycle.
  
  set RTOL [expr {$R + $TOL}]

  # restore the complete list of atoms to consider:
  
  if { $step % 100 == 0 } { cleardrops }	

  while {[nextatom]} { 
  
    set rvec [getcoord]
    foreach { x y z } $rvec { break }
    set rho [expr {sqrt(($x-$x0)*($x-$x0) + ($y-$y0)*($y-$y0) + \
      ($z-$z0)*($z-$z0))}]
	  
    # Atoms at distances 0 to $R from the bubble center are pushed,
    # atoms father than $RTOL are dropped. Atoms between $R to $RTOL,
    # that is with a layer of thickness $TOL, are neither pushed nor 
    # dropped at this step, so that they will be considered again at 
    # the next step(s). They may come closer to the bubble and then 
    # they will have to be pushed.

    if { $rho < $R } {
	
      set forceX [expr {$K * ($x-$x0) / $rho}]
      set forceY [expr {$K * ($y-$y0) / $rho}]
      set forceZ [expr {$K * ($z-$z0) / $rho}]
      addforce "$forceX $forceY $forceZ"

    } elseif { $rho > $RTOL } {
	
	  dropatom

    }
  }
}
