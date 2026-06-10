# --------------------------------------------------------------------
# HELPER FUNCTION: Find files with a prefix, suffix and an optional number within a given range
# Arguments:
#   dir_path			 :  directory path to search
#   prefix and suffix    :  file name prefix and suffix
#   min_num and max_num  :  optional range (both inclusive). "" for none

proc find_files {dir_path prefix suffix {min_num ""} {max_num ""} {sort_natural 1} {return_abs_path 0}} {
    #if {![file isdirectory $dir_path]} { error "Directory '$dir_path' not found." }
    set result_list {}; set pattern "${prefix}(\[0-9\]+)${suffix}$";
    foreach f [glob -nocomplain -directory $dir_path *] {
        if {[file isfile $f]} { set filename [file tail $f];
            if {[regexp $pattern $filename -> num]} { set num [expr {$num + 0}];    # ensure numeric
                if {($min_num eq "" || $num >= $min_num) && ($max_num eq "" || $num <= $max_num)} {
                    if { $return_abs_path == 1 } { set fpath [file normalize $f]; } else { set fpath [file join $dir_path $filename]; }
					lappend result_list $fpath; }}}}
    if { $sort_natural == 1 } { set result_list [lsort -dictionary $result_list] };
    return $result_list;
}

# Read atomselection string from file (must contain only a single line)
proc read_selection_from_file { file_path } {
	set fd [open $file_path];
	set data [read $fd];
	close $fd;

	set trimmed [string trim $data];
	return $trimmed;
}
