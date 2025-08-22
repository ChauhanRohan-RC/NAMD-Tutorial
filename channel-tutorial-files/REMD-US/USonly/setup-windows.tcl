set num_replicas 19

for {set i 0} {$i < $num_replicas} {incr i} {

## write the config files
  set in1 [open "win-base.conf" r]
  set out [open "example-output/${i}/win${i}A.conf" w]
  puts $out "set num $i"
  while { [gets $in1 line] >= 0} {
    puts $out $line
  } 
  close $out

## write the colvars input files      
  set in2 [open "US-base.in" r]
  set out [open "example-output/${i}/US-win${i}.in" w]
  while { [gets $in2 line] >= 0} {
    if { [string match "*CENTER*" $line]} {
      puts $out "   centers        [expr 5 - $i]"
    } else {  
      puts $out $line
    }
  }
}

