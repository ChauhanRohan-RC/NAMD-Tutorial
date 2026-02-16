#!/usr/bin/env -S vmd -dispdev text -e

###########################################################################################
# Script to calculate RMSF of each residue of a molecule in VMD						 ##
###########################################################################################
# NOTE: plotting data requires "plot_data.py"

# **REQUIRE bigdcd.tcl**
package require bigdcd;
# OR
# source bigdcd.tcl;

### USAGE-------------------
# 1. Copy this script in working dir
# 2. INPUT: Set input structure (.psf) and frame (.dcd, .pdb, .coor) files
# 3. INPUT: Set BASE_SELECTION and RESIDUE_SELECTION
# 4. run with "./rmsf.tcl"
# 5. Generates output file "rmsf.csv" and "rmsf.pdf" (if plot_output = 1)

## NOTE UNITS: RMSF values are in Angstrom (Å)

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
# --------------------------------------------------------------------


set COL_NAME_RESIDUE	"RESIDUE";	# Frame index
set COL_NAME_RMSF		"RMSF";		# RMSF (Å)


# ==================
# INPUT
# ==================
set psf_file		"../../common/amyld_wb.psf";		# TODO: input strcuture file (.psf)

# TODO: LIST of trajectory (.dcd) or single frame (.pdb, .coor) files separated by space	
# set frame_files		{ "../dna_wb_eq.dcd" "../dna_wb_eq1.dcd" };
set frame_files	[find_files ".." "amyld_wb_eq" ".dcd"];  # <dir> <prefix> <suffix> [min_num] [max_num]

set frame_index_ref		0;		# Reference Frame Index (default 0)
set frame_index_start 	-1;		# Inclusive, -1 for None
set frame_index_end 	-1;		# Exclusive, -1 for None


# TODO: BASE SELECTION for each Frame.
# -> Used to select the part of system conating residues
# -> Examples: "protein", "nucleic", "lipid", "all not water"
set base_selection		"protein";

# TODO: RESIDUE SELECTION within BASE_SELECTION
# -> Used to get RESIDUE_ID's of residues within BASE_SELECTION for which RMSF is required
# -> EXAMPLES
#    1. for Protein all residues: "protein and alpha" or "protein and name CA" 
#    2. for Nucleic all residues: "nucleic and name C1'"
#    3. for RESIDUE RANGE, append: "and (resid <start> to <end>)" or "and (resid <one> <two> <three> ...)"
#       Example: "protein and alpha and (resid 10 to 14)", "nucleic and name C1' and not (resid 11 15)"
set residue_selection 	"${base_selection} and name CA"


# ==================
# INPUT
# ==================
set out_file_prefix		"rmsf"
set out_data_file 		"${out_file_prefix}.csv";		# output data file

# Format
set out_format			"%.6f";		# RG Format.
set out_delimiter 		" \t ";		# output delimiter

# Plotting
set plot_output				1;				# [0/1] plot output using "plot_data.py"
set show_interactive_plot 	1; 				# [0/1]
set out_plot_file 			"${out_file_prefix}.pdf";	# [Optional] Only if plot_output = 1

set comment_token 		"#";	# Token used for Comments
set comment_header		0;		# Whether to comment out the columns header
# -----------------------------------------------------



# ==============================================
# MAIN
# ==============================================

puts "\n=============================================="
puts "=============  RMSF (bigdcd)  ================="
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
if { [info exists frame_index_ref] == 0 || $frame_index_ref < 0 } {
	puts "\n---------------------------------------------------------"
	puts " => ERROR: Reference Frame Index must be defined and >= 0"
	puts "---------------------------------------------------------\n"
	exit;
}

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

set sel_cur [atomselect $mol_id "$base_selection"];
set num_atoms_sel [$sel_cur num];
if { $num_atoms_sel == 0 } {
	puts "\n-------------------------------------------------------------"
	puts " => ERROR: BASE_SELECTION has no atoms : \"$base_selection\""
	puts "-------------------------------------------------------------\n"
	exit;
}


set res_ids [[atomselect $mol_id "$residue_selection"] get resid];
set res_count [llength $res_ids];

if { $res_count == 0 } {
	puts "\n-------------------------------------------------------------"
	puts " => ERROR: RESIDUE_SELECTION has no residues : \"$residue_selection\""
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
	global found_ref_frame ref_frame_filename;

	if { [info exists out_file] == 1 } {
		flush $out_file;
		close $out_file;
	}

	if { [info exists mol_id] == 1 } {
		mol delete $mol_id;
	}

	if { [info exists found_ref_frame] == 1 && $found_ref_frame == 1} {
		file delete $ref_frame_filename;
	}
}


# Functions ---------------------------------

# Reference Frame vars
set ref_frame_filename "rmsd_ref_frame-[clock seconds].temp.coor";
set sel_ref 	"";
set found_ref_frame 0;
set num_frames 0;

# Bigdcd dryrun to find and save ref frame
proc try_save_ref_frame_bigdcd { i } {
	global frame_index_ref;

	set i [expr $i - 1];

	if { $i == $frame_index_ref } {
		global mol_id allatoms;
		global found_ref_frame ref_frame_filename;

		# Saving ref frame to pdb file
# 		set alla [atomselect $mol_id "all" frame $i];
		$allatoms writenamdbin $ref_frame_filename;
		set found_ref_frame 1;

		puts "\n-------------------------------------------------------------"
		puts " => bigdcd DRY-RUN: Reference Frame Found at Index : ${i}"
		puts "-------------------------------------------------------------\n"
		error "--";		# Exit bigdcd dryrun
	}
}


# Main function to calc rmsf using bigdcd
proc calc_rmsf { i } {
	global frame_index_ref frame_index_start frame_index_end num_frames;
	global base_selection
	global mol_id sel_ref sel_cur allatoms;
	global res_ids rmsf_arr sel_ref_res_arr sel_cur_res_arr;

	set i [expr $i - 1];	# required for bigdcd

	if { $i < $frame_index_start } {
		return;
	}

	if { $frame_index_end > 0 && $i >= $frame_index_end } {
		# Finished. Exit bigdcd
		error "Finished: FRAME RANGE \[$frame_index_start, $frame_index_end) processed";
	}

	if {[expr $i % 10000] == 0} {
		puts "\n------------------------------"
		puts " => processing Frame ${i}";
		puts "------------------------------\n"
	}

	## Note: bigdcd already does this
	#$allatoms frame $i;
	#$sel_cur frame $i;
	#$sel_ref frame 0;

	# Alignment
	$sel_cur move [measure fit $sel_cur $sel_ref];

	# compute the contribution to RMSF (loop through all residues)
	foreach r $res_ids {

		set ref $sel_ref_res_arr($r);
		set cur $sel_cur_res_arr($r);

		set rmsd_r [measure rmsd $cur $ref]
		set rmsf_arr($r) [ expr $rmsf_arr($r) + ($rmsd_r * $rmsd_r) ]
	}

	incr num_frames;
}


# =====================================
# DRY-RUN: Find Ref Frame
# =====================================
eval "bigdcd try_save_ref_frame_bigdcd auto [join $frame_files]";
bigdcd_wait;

if { $found_ref_frame == 0 } {
	puts "\n-------------------------------------------------------------"
	puts " => ERROR: Could not find Reference Frame Index : $frame_index_ref"
	puts "-------------------------------------------------------------\n"
	cleanup;
	exit;
}

# Setup fresh molecule for main run
mol delete $mol_id;

set mol_id [mol new $psf_file];
mol addfile "$ref_frame_filename" type namdbin waitfor all molid $mol_id;
set sel_ref [atomselect $mol_id "$base_selection" frame 0];
set sel_cur [atomselect $mol_id "$base_selection"];
set allatoms [atomselect $mol_id "all"];

# Initialize Output Arrays
foreach r $res_ids {
	set rmsf_arr($r) 0;	# Output array

	set sel_str "${base_selection} and resid $r and noh";
	set sel_ref_res_arr($r) [atomselect $mol_id $sel_str frame 0];
	set sel_cur_res_arr($r) [atomselect $mol_id $sel_str];
}


# =====================================
# MAIN-RUN: Calculate RMSF
# =====================================
eval "bigdcd calc_rmsf auto [join $frame_files]";
bigdcd_wait;

if { $num_frames ==0 } {
	puts "\n-------------------------------------------------------------"
	puts " => ERROR: No Frame Processed"
	puts "-------------------------------------------------------------\n"
	cleanup;
	exit;
}

# Final RMSF Calculation ------------
set average_rmsf 0;

foreach r $res_ids {
	# Final RMSF Calculation
    set rmsf_arr($r) [ expr sqrt($rmsf_arr($r) / $num_frames) ]

    # Computing average RMSF
    set average_rmsf [expr $average_rmsf + $rmsf_arr($r)]

    # Set the "user" field (of atoms of the residue) to RMSF value
#     set res_atoms $sel_ref_res_arr($r)
#     $res_atoms set user $rmsf_arr($r);
}

set average_rmsf [ expr $average_rmsf / $res_count ]



# ===============================
# OUTPUT File
# ===============================
set out_file [open $out_data_file w]

log "${comment_token} ==================== RMSF (bigdcd) ======================"
log "${comment_token} INPUT Structure File: \"${psf_file}\" | Frame File(s): \[${frame_files}\]"
log "${comment_token} INPUT Base Selection: \"${base_selection}\""
log "${comment_token} INPUT Residue Selection: \"${residue_selection}\""
log "${comment_token} => ATOM COUNT: Total ($num_atoms_all), Base Selection ($num_atoms_sel)"
log "${comment_token} => Selected RESIDUE COUNT: $res_count"
log "${comment_token}-----------------------------------------------------------"
log "${comment_token} INPUT REFERENCE Frame Index: ${frame_index_ref}";
log "${comment_token} INPUT FRAME Index RANGE: \[${frame_index_start}, ${frame_end_str})"
log "${comment_token} => Frames Processed: ${num_frames}";
log "${comment_token}-----------------------------------------------------------"
log "${comment_token} AVERAGE RMSF: [format ${out_format} ${average_rmsf}] Å"
log "${comment_token} Units: RMSF (Angstrom Å)"
log "${comment_token}-----------------------------------------------------------"

# Header
set out_header "${COL_NAME_RESIDUE}${out_delimiter}${COL_NAME_RMSF}";
if {$comment_header == 1} {
	set out_header "${comment_token}${out_header}";
}

log2file $out_header;

# RMSF DATA
foreach r $res_ids {
	log "${r}${out_delimiter}[format ${out_format} $rmsf_arr($r)]"
}


# GC ----------
cleanup;

# Plot Output -------------------
if { $plot_output == 1 } {
	set cmd "./plot_data.py -t {RMSF Plot} -xl {Residue} -yl {RMSF (Å)} -o \"${out_plot_file}\" \"${out_data_file}\"";
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



