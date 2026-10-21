#!/usr/bin/env -S vmd -dispdev text -e

# =========================================================================================
# Script to find Buried Surface Area (BSA) of a dimer
# => BSA = SASA (chain_1) + SASA (chain_1) - SASA (dimer)
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
set psf_file		"../../common/amyld_wb.psf";					# TODO: input strcuture file (.psf)

# TODO: Trajectory Files (.pdb, .coor, .dcd etc)
set frame_files	    [find_files "../" "amyld_wb_eq" ".dcd"];  # <dir> <prefix> <suffix> [min_num] [max_num]

set frame_index_start 	-1;		# Inclusive, -1 for None
set frame_index_end 	-1;		# Exclusive, -1 for None

set selection_dimer			"protein"
set selection_chain1		"protein and chain A";		# TODO: 
set selection_chain2		"protein and chain B";		# TODO: 

set sasa_probe_radius		1.4;	# 1.4 is standard for water

# -----------------------
# OUTPUT
# -----------------------
set out_data_file		"buried_sasa.csv"

set out_sep 			" "
set comment_token 		"#";	# Token used for Comments
set comment_header		0;		# Whether to comment out the columns header

set progress_log_frames 1000;


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

set sel_dimer [atomselect $mol_id "$selection_dimer"]
set sel_chain1 [atomselect $mol_id "$selection_chain1"];
set sel_chain2 [atomselect $mol_id "$selection_chain2"];
set num_atoms_dimer [$sel_dimer num];
set num_atoms_chain1 [$sel_chain1 num];
set num_atoms_chain2 [$sel_chain2 num];

if { $num_atoms_chain1 == 0 } {
	puts "\n-------------------------------------------------------------"
	puts " => ERROR: Chain 1 has no atoms : \"$selection_chain1\""
	puts "-------------------------------------------------------------\n"
	exit;
}

if { $num_atoms_chain2 == 0 } {
	puts "\n-------------------------------------------------------------"
	puts " => ERROR: Chain 2 has no atoms : \"$selection_chain2\""
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


set out_file [open $out_data_file w];

log "${comment_token} ========= Buried Surface Area (bigdcd) ==========="
log "${comment_token} INPUT Structure File      : \"${psf_file}\""
log "${comment_token} INPUT Frame File(s)       : \[${frame_files}\]"
log "${comment_token}-----------------------------------------------------------"
log "${comment_token} INPUT SASA Probe Radius   : ${sasa_probe_radius} Å"
log "${comment_token} INPUT Selection Dimer     : \"${selection_dimer}\""
log "${comment_token} INPUT Selection Chain 1   : \"${selection_chain1}\""
log "${comment_token} INPUT Selection Chain 2   : \"${selection_chain2}\""
log "${comment_token} => TOTAL Atoms  			: $num_atoms_all"
log "${comment_token} => TOTAL Dimer Atoms      : $num_atoms_dimer"
log "${comment_token}    => Chain 1 Atoms       : $num_atoms_chain1"
log "${comment_token}    => Chain 2 Atoms       : $num_atoms_chain2"
log "${comment_token}-----------------------------------------------------------"
log "${comment_token} INPUT FRAME Index RANGE: \[${frame_index_start}, ${frame_end_str})";
#log "${comment_token} => Frames Processed: ${num_frames}";		# num_frames not initialized yet
log "${comment_token}-----------------------------------------------------------";
log "${comment_token} BSA (buried Surface area): (SASA_CHAIN1 + SASA_CHAIN2) - SASA_PROT"
log "${comment_token} UNITS: Surface Area => Å^2"
log "${comment_token}-----------------------------------------------------------";

# Header
set out_header "FRAME${out_sep}SASA_CHAIN1${out_sep}SASA_CHAIN2${out_sep}SASA_PROT${out_sep}BSA";
if {$comment_header == 1} {
	set out_header "${comment_token}${out_header}";
}

log2file $out_header;


# Functions ---------------------------------

set num_frames 0;


# Main function to calc SASA
proc calc_sasa { i } {
	global frame_index_start frame_index_end num_frames;
	global mol_id sel_dimer sel_chain1 sel_chain2;
	global out_sep progress_log_frames;
	global sasa_probe_radius;
	
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

    ## Note: bigdcd already does this
    #$sel_chain1 frame $i
    #$sel_chain2 frame $i
    #$sel update

    set sasa_dimer [measure sasa $sasa_probe_radius $sel_dimer -restrict $sel_dimer];
    set sasa_chain1 [measure sasa $sasa_probe_radius $sel_chain1 -restrict $sel_chain1];
    set sasa_chain2 [measure sasa $sasa_probe_radius $sel_chain2 -restrict $sel_chain2];
    set bsa [expr ($sasa_chain1 + $sasa_chain2) - $sasa_dimer]

    # Format the numbers to 2 decimal places and write the row to the file
    set out_line [format "%d${out_sep}%.2f${out_sep}%.2f${out_sep}%.2f${out_sep}%.2f" $i $sasa_chain1 $sasa_chain2 $sasa_dimer $bsa]
    log2file $out_line
    
    incr num_frames;
}




# ================================================
# MAIN-RUN: Calculate Sec Structure Percentage
# ================================================
eval "bigdcd calc_sasa auto [join $frame_files]";
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

puts "\n================= BSA FINISHED  ==================="
puts "=> FRAMES processed : $num_frames"
puts "-> Time Taken       : [expr $time_end -$time_start] secs"
puts "=> OUTPUT Data File : \"${out_data_file}\""
puts "===================================================\n"

exit;

