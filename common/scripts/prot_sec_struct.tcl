#!/usr/bin/env -S vmd -dispdev text -e

# =========================================================================================
# Finds the secondary structure percentages of Protein over trajectory file(s)
# =========================================================================================

# **REQUIRE bigdcd.tcl**
package require bigdcd;
# OR
# source bigdcd.tcl;

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




# -----------------------
# INPUT
# -----------------------
set psf_file		"../../common/amyl_wb.psf";					# TODO: input strcuture file (.psf)

# TODO: Trajectory Files (.pdb, .coor, .dcd etc)
set frame_files	    [find_files "../" "amyl_wb_eq" ".dcd"];  # <dir> <prefix> <suffix> [min_num] [max_num]

set frame_index_start 	-1;		# Inclusive, -1 for None
set frame_index_end 	-1;		# Exclusive, -1 for None

set residue_selection		"protein and name CA";		# TODO: protein residue selection (1 alpha-C per residue)

# -----------------------
# OUTPUT
# -----------------------
set out_structure_file		"prot_sec_structure.txt"

set out_sep 			" "
set comment_token 		"#";	# Token used for Comments
set comment_header		0;		# Whether to comment out the columns header

set progress_log_frames 100;


# -----------------------
# MAIN
# -----------------------

puts "\n=========================================================="
puts "=========  Protein Secondary Structure (bigdcd)  ========="
puts "==========================================================\n"

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

# Selection -----------
set num_atoms_all [molinfo $mol_id get numatoms];
if { $num_atoms_all == 0 } {
	puts "\n-------------------------------------------------------------"
	puts " => ERROR: Molecule has no atoms : \"$psf_file\""
	puts "-------------------------------------------------------------\n"
	exit;
}

set sel_cur [atomselect $mol_id "$residue_selection"];
set res_count [$sel_cur num];
if { $res_count == 0 } {
	puts "\n-------------------------------------------------------------"
	puts " => ERROR: No residues found in selection : \"$residue_selection\""
	puts "-------------------------------------------------------------\n"
	exit;
}


# =================================
# Output File
# =================================

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


set out_file [open $out_structure_file w];

log "${comment_token} ========= Protein Secondary Structure (bigdcd) ==========="
log "${comment_token} INPUT Structure File      : \"${psf_file}\""
log "${comment_token} INPUT Frame File(s)       : \[${frame_files}\]"
log "${comment_token}-----------------------------------------------------------"
log "${comment_token} INPUT Residue Selection   : \"${residue_selection}\""
log "${comment_token} => TOTAL Atoms  			: $num_atoms_all"
log "${comment_token} => Selected RESIDUE COUNT : $res_count"
log "${comment_token}-----------------------------------------------------------"
log "${comment_token} INPUT FRAME Index RANGE: \[${frame_index_start}, ${frame_end_str})";
#log "${comment_token} => Frames Processed: ${num_frames}";		# num_frames not initialized yet
log "${comment_token}-----------------------------------------------------------";
log "${comment_token} UNITS: Percentage (%)"
log "${comment_token}-----------------------------------------------------------";

# Header
set out_header "Frame${out_sep}HELIX${out_sep}SHEET${out_sep}TURN${out_sep}COIL${out_sep}OTHER";
if {$comment_header == 1} {
	set out_header "${comment_token}${out_header}";
}

log2file $out_header;


# Functions ---------------------------------

set num_frames 0;


# Main function to calc sec_struct using bigdcd
proc calc_ss { i } {
	global frame_index_start frame_index_end num_frames;
	global mol_id sel_cur res_count;
	global out_sep progress_log_frames;
	
	set i [expr $i - 1];	# required for bigdcd

	if { $i < $frame_index_start } {
		return;
	}

	if { $frame_index_end > 0 && $i >= $frame_index_end } {
		# Finished. Exit bigdcd
		error "Finished: FRAME RANGE \[$frame_index_start, $frame_index_end) processed";
	}

	if {[expr $i % $progress_log_frames] == 0} {
		puts "\n------------------------------"
		puts " => processing Frame ${i}";
		puts "------------------------------\n"
	}

    ## Go to the current frame and update the display
    #animate goto $i
    #display update

    ## Force VMD to recalculate secondary structure for this specific frame via STRIDE
    vmd_calculate_structure $mol_id

    ## Note: bigdcd already does this
    #$sel_cur frame $i
    #$sel update

    # Get the secondary structure code for every residue
    set ss_list [$sel_cur get structure]

    # Initialize counters
    set count_helix 0
    set count_sheet 0
    set count_turn  0
    set count_coil  0
    set count_other 0

    ## Count the occurrences of each structure type
    # STRIDE Codes:
    # H (Alpha Helix), G (3_10 Helix), I (Pi Helix) -> Helix
    # E (Extended Sheet), B (Isolated Bridge) -> Sheet
    # T -> Turn
    # C -> Coil
    foreach ss $ss_list {
        if {$ss == "H" || $ss == "G" || $ss == "I"} {
            incr count_helix
        } elseif {$ss == "E" || $ss == "B"} {
            incr count_sheet
        } elseif {$ss == "T"} {
            incr count_turn
        } elseif {$ss == "C"} {
            incr count_coil
        } else {
            incr count_other
        }
    }

    # Calculate percentages (multiplied by 100.0 to enforce floating-point division)
    set p_helix [expr {($count_helix * 100.0) / $res_count}]
    set p_sheet [expr {($count_sheet * 100.0) / $res_count}]
    set p_turn  [expr {($count_turn * 100.0) / $res_count}]
    set p_coil  [expr {($count_coil * 100.0) / $res_count}]
    set p_other [expr {($count_other * 100.0) / $res_count}]

    # Format the numbers to 2 decimal places and write the row to the file
    set out_line [format "%d${out_sep}%.2f${out_sep}%.2f${out_sep}%.2f${out_sep}%.2f${out_sep}%.2f" $i $p_helix $p_sheet $p_turn $p_coil $p_other]
    log2file $out_line
    
    incr num_frames;
}




# ================================================
# MAIN-RUN: Calculate Sec Structure Percentage
# ================================================
eval "bigdcd calc_ss auto [join $frame_files]";
bigdcd_wait;

# GC ----------
cleanup;

if { $num_frames == 0 } {
	puts "\n-------------------------------------------------------------"
	puts " => ERROR: No Frame Processed"
	puts "-------------------------------------------------------------\n"
	exit;
}


set time_end [clock seconds];

puts "\n========= Protein Sec Structure FINISHED  ========="
puts "=> OUTPUT Data File : \"${out_structure_file}\""
puts "-> Time Taken       : [expr $time_end -$time_start] secs"
puts "===================================================\n"

exit;

