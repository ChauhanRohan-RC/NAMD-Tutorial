#!/usr/bin/env -S vmd -dispdev text -e

###############################################################################
# Script to calculate Dihedral angles from trajectory in VMD         		 ##
###############################################################################
# NOTE: Plotting data requires "plot_data.py"

# **REQUIRE bigdcd.tcl**
package require bigdcd;
# OR
# source bigdcd.tcl;

### USAGE ------------------
# 1. Copy this script in working dir
# 2. INPUT: Set input structure (.psf) and frame (.dcd, .pdb, .coor) files
# 3. INPUT: Set ATOM_SELECTION
# 4. run with "./dihedral.tcl"
# 5. Generates output file "dihedral.csv" and "dihedral.pdf" (if plot_output = 1)

## NOTE UNITS: Dihedral Angle values are in degrees (°)


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
# --------------------------------------------------------------------


set COL_NAME_FRAME		"FRAME";		# Frame index
set COL_NAME_DIHEDRAL	"DIHEDRAL";		# Dihedral Angle (degrees)


# ======================
# INPUT
# ==================
set psf_file		"../../common/amyld_wb.psf";		# TODO: input strcuture file (.psf)

# TODO: LIST of trajectory (.dcd) or single frame (.pdb, .coor) files separated by space
#set frame_files		{ "../dna_wb_eq.dcd" "../dna_wb_eq1.dcd" };
set frame_files	[find_files ".." "amyld_wb_eq" ".dcd"];  # <dir> <prefix> <suffix> [min_num] [max_num]

set frame_index_start 	-1;		# Inclusive, -1 for None
set frame_index_end 	-1;		# Exclusive, -1 for None

# TODO: ATOM SELECTION's for dihedral angle
set atom_sel1		"protein and resid 1 and name CA";
set atom_sel2		"protein and resid 1 and name CB";
set atom_sel3		"protein and resid 1 and name CG";
set atom_sel4		"protein and resid 1 and name OD1";


# ==================
# OUTPUT
# ==================
set out_file_prefix		"dihedral"
set out_data_file 		"${out_file_prefix}.csv";		# output data file

# Format
set out_format			"%.6f";		# dihedral Format.
set out_delimiter 		" ";		# output delimiter

# Plotting
set plot_output				1;				# [0/1] plot output using "plot_data.py"
set show_interactive_plot 	1; 				# [0/1]
set out_plot_file 			"${out_file_prefix}.pdf";	# [Optional] Only if plot_output = 1


set comment_token 		"#";	# Token used for Comments
set comment_header		0;		# Whether to comment out the columns header



# ==============================================
# MAIN
# ==============================================

puts "\n=============================================="
puts "========  Dihedral Angle (bigdcd)  ============"
puts "==============================================\n"

set time_start [clock seconds];

## CHECKS ==========================

# Frame files -----------
if {[info exists frame_files] == 0 || [llength $frame_files] == 0} {
	puts "\n------------------------------------------"
	puts " => ERROR: No Frame Files specified !!"
	puts "------------------------------------------\n"
	exit;
}


# Frame Range ---------------------
if { [info exists frame_index_start] == 0 || $frame_index_start < 0 } {
	set frame_index_start 0;
}

set frame_end_str ""
set frame_count_str ""
if { [info exists frame_index_end] == 0 || $frame_index_end < 0 } {
	set frame_index_end -1;
	set frame_end_str "LAST"
	if { $frame_index_start > 0 } {
		set frame_count_str "$frame_index_start-LAST"
	} else {
		set frame_count_str "ALL"
	}
} else {
	if { $frame_index_end <= $frame_index_start } {
		puts "\n------------------------------------------"
		puts " => ERROR: FRAME_INDEX_END must be > FRAME_INDEX_START. given start: $frame_index_start, end: $frame_index_end"
		puts "------------------------------------------\n"
		exit;
	}

	set frame_end_str "$frame_index_end"
	set frame_count_str "[expr $frame_index_end - $frame_index_start]"
}



# Loading Molecule -------------
set mol_id [mol new $psf_file waitfor all];

# foreach frame_file $frame_files {
# 	puts "---------------------------------------------"
# 	puts " -> Loading Frames from file : ${frame_file}"
# 	puts "---------------------------------------------"
# 	mol addfile $frame_file waitfor all molid $mol_id;
# }

# Selection -----------
set allatoms [atomselect $mol_id "all"];
set num_atoms_all [$allatoms num];
if { $num_atoms_all == 0 } {
	puts "\n-------------------------------------------------------------"
	puts " => ERROR: Molecule has no atoms : \"$psf_file\""
	puts "-------------------------------------------------------------\n"
	exit;
}


proc validate_sel { sel } {
	# No atoms
	if { [$sel num] == 0 } {
		puts "\n-------------------------------------------------------------"
		puts " => ERROR: Selection has no atoms : \"[$sel text]\""
		puts "-------------------------------------------------------------\n"
		return 0;
	}

	# More than 1 atom
	if { [$sel num] > 1 } {
		puts "\n-------------------------------------------------------------"
		puts " => ERROR: Selection has more than 1 atom : \"[$sel text]\" (atom count: [$sel num])"
		puts "-------------------------------------------------------------\n"
		return 0;
	}

	return 1;
}


# Atom Selections
set sel1 [atomselect $mol_id "$atom_sel1"];
set sel2 [atomselect $mol_id "$atom_sel2"];
set sel3 [atomselect $mol_id "$atom_sel3"];
set sel4 [atomselect $mol_id "$atom_sel4"];

# Validation
set val_failed 0;
set atom_indices [list];

foreach sel [list $sel1 $sel2 $sel3 $sel4] {
	if { [validate_sel $sel] == 0 } {
		set val_failed 1;
# 		break;	# check all selections regardless
	}

	lappend atom_indices [lindex [$sel get index] 0];
}

if { $val_failed == 1 } {
	exit;
}





# =================================
# Output File
# =================================

set out_file [open $out_data_file w];

proc log2file { msg } {
	global out_file;
	puts $out_file $msg;	# to output file
}

proc log { msg } {
	global out_file;
	puts $msg; flush stdout;		# to stdout
	log2file $msg;					# to output file
}

proc cleanup { } {
	global out_file mol_id;

	if { [info exists out_file] == 1 } {
		flush $out_file;
		close $out_file;
	}

	if { [info exists mol_id] == 1 } {
		mol delete $mol_id;
	}
}

log2file "${comment_token} ================ Dihedral Angle (bigdcd) =================="
puts "----------------------------------------------"
log "${comment_token} INPUT Structure File: \"${psf_file}\" | Frame File(s): \[${frame_files}\]"
log "${comment_token} INPUT Atom Selection-1: \"${atom_sel1}\""
log "${comment_token} INPUT Atom Selection-2: \"${atom_sel2}\""
log "${comment_token} INPUT Atom Selection-3: \"${atom_sel3}\""
log "${comment_token} INPUT Atom Selection-4: \"${atom_sel4}\""
log "${comment_token} => ATOM COUNT: Total ($num_atoms_all)"
log "${comment_token}-----------------------------------------------------------"
log "${comment_token} FRAME Index RANGE: \[${frame_index_start}, ${frame_end_str})  |  COUNT: $frame_count_str"
log "${comment_token} Units: Dihedral Angle (degrees °, in range -180° to +180°)"
log "${comment_token}-----------------------------------------------------------"

set out_header "${COL_NAME_FRAME}${out_delimiter}${COL_NAME_DIHEDRAL}";
if {$comment_header == 1} {
	set out_header "${comment_token}${out_header}";
}

log2file $out_header;
#puts $out_header;


# =====================================
# MAIN-RUN: Calculate Dihedral
# =====================================

proc calc_dihed { i } {
	global frame_index_start frame_index_end;
	global mol_id atom_indices allatoms;
	global out_delimiter out_format;

	set i [expr $i - 1];

	if { $i < $frame_index_start } {
		return;
	}

	if { $frame_index_end > 0 && $i >= $frame_index_end } {
		# Finished. Exit bigdcd
		error "Finished: FRAME RANGE \[$frame_index_start, $frame_index_end) processed";
	}


# 	$sel_cur frame $i;	# done by bigdcd
	set dihed_val [measure dihed $atom_indices molid $mol_id];
	set dihed_str [format "$out_format" $dihed_val];

	if {[expr $i % 10000] == 0} {
		puts "\n-----------------------------------------------"
		puts " => processing Frame $i -> DIHED: ${dihed_str}°";
		puts "-----------------------------------------------\n"
	}

	set out_line "${i}${out_delimiter}${dihed_str}";
	log2file $out_line;
}

# Execute Calculation (bigdcd)
eval "bigdcd calc_dihed auto [join $frame_files]";
bigdcd_wait;


# GC
cleanup;

# Plot Output -------------------
if { $plot_output == 1 } {
	set cmd "./plot_data.py -t {Dihedral Angle} -xl {Frame} -yl {DIHED (°)} -o \"${out_plot_file}\" \"${out_data_file}\"";
	if { $show_interactive_plot == 0 } {
		set cmd "$cmd -ni";
	}

	eval $cmd;
}


set time_end [clock seconds];

puts "\n==================  FINISHED  ====================="
puts "=> OUTPUT Data File: \"${out_data_file}\""
if { $plot_output == 1 } {
	puts "=> OUTPUT Plot File: \"${out_plot_file}\""
}
puts "-> Time Taken: [expr $time_end -$time_start] secs"
puts "===================================================\n"

exit;
