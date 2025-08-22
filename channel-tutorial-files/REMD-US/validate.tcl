set num_replicas 5

proc replica_neighbors { i } {
  set n 5
  if { $i % 2 } { set s -1 } { set s 1 }
  set result {}
  foreach dx { $s 0 -$s } {
    set	j [expr	$i + $dx]
    if { $j < 0 || $j >= $n } {
      lappend result $i ;# swap	with self
    } {
      lappend result $j	;# swap	with neighbor
    }
  }
  return $result
}

for { set i 0 } { $i < $num_replicas } { incr i } {
  set j 0
  foreach nbr [replica_neighbors $i] {
    if { $nbr < 0 } {
      error "replica_neighbors inconsistency detected: neighbor $j of replica $i is $nbr but should not be negative"
    }
    if { $nbr >= $num_replicas } {
      error "replica_neighbors inconsistency detected: neighbor $j of replica $i is $nbr but there are only $num_replicas replicas"
    }
    set rnbrl [replica_neighbors $nbr]
    set rnbrc [llength $rnbrl]
    if { $j >= $rnbrc } {
      error "replica_neighbors inconsistency detected: neighbor $j of replica $i is $nbr but replica $nbr has only $rnbrc neighbors"
    }
    set rnbr [lindex $rnbrl $j]
    if { $rnbr != $i } {
      error "replica_neighbors inconsistency detected: neighbor $j of replica $i is $nbr but neighbor $j of replica $nbr is $rnbr"
    }
    incr j
  }
}
