# sphere.tcl -- make a linear potential leading to a flat spherical well

### PARAMETERS ###
#
set gridlen 30.0
set gridnum 30
set center "0.0 0.0 10.0"
set radius 5.0
set gradient 1.0
set outfile potassium.dx
#set outfile chloride.dx


### MAIN ###
#

# calculate number of gridpoints
set gridorigin [vecadd [vecscale "$gridlen $gridlen $gridlen" -0.5] $center]
set gridpoints [expr {$gridnum * $gridnum * $gridnum}]
set griddelta [expr {$gridlen/$gridnum}]

# open output file
if { [file exists $outfile] } {
    puts "File $outfile already exists!"
    #exit
}
set out [open $outfile w]

# write DX header
puts $out "object 1 class gridpositions counts $gridnum $gridnum $gridnum"
puts $out "origin $gridorigin"
puts $out "delta $griddelta 0.0 0.0"
puts $out "delta 0.0 $griddelta 0.0"
puts $out "delta 0.0 0.0 $griddelta"
puts $out "object 2 class gridconnections counts $gridnum $gridnum $gridnum"
puts $out "object 3 class array type double rank 0 items $gridpoints data follows"

# calculate and write grid values
set n 0		;# counter to print 3 values per line
for { set i 0 } { $i < $gridlen } { incr i } {
    set x [expr {$i * $griddelta + [lindex $gridorigin 0]}]
    for { set j 0 } { $j < $gridlen } { incr j } {
	set y [expr {$j * $griddelta + [lindex $gridorigin 1]}]
	for { set k 0 } { $k < $gridlen } { incr k } {
	    set z [expr {$k * $griddelta + [lindex $gridorigin 2]}]
	    set r [veclength [vecsub "$x $y $z" $center]]
	    if { $r > $radius } {
		set v [expr {($r - $radius) * $gradient}]
	    } else {
		set v 0
	    }
	    puts -nonewline $out $v
	    incr n
	    if { $n == 3 } {
		# print newline and reset counter
		puts $out ""
		set n 0
	    } else {
		puts -nonewline $out " "
	    }
	}
    }
}

close $out

exit