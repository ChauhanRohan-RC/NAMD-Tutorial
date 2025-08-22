set num_replicas 19
set ns_per_job 0.5
set t [lindex $argv 0]
set eqt 1.0
set ni [expr int($eqt / $ns_per_job)]
set nf [expr int($t / $ns_per_job)]

set INPUT [open "INPUT_1-5ns" w]

puts $INPUT "\#\# wham -13 5 180 0.0001 310.0 0 INPUT_1-5ns AmtB_REMD_US_1-5ns.pmf"
puts $INPUT "\#\#/path/to/timeseries/file loc_win_min spring \[correl time\] \[temp\]"
puts $INPUT "\#\#"

for {set i 0} {$i < $num_replicas} {incr i} {
	for { set j $ni } { $j < $nf } { incr j } {
		puts $INPUT "$i/AmtB-A.job$j.$i.sort.100.colvars.traj [expr 5.0 - $i] 2.5"
	}
}

close $INPUT
