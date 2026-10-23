#!/usr/bin/env -S vmd -dispdev text -e

#############################################################################
# Script to calculate RMSD, RMSF and RG from trajectory files in VMD	   ##
#############################################################################
# NOTE: Plotting data requires "plot_data.py"

# **REQUIRE bigdcd.tcl**
package require bigdcd;
# OR
# source bigdcd.tcl;

### USAGE ------------------
# 1. Copy this script in working dir
# 2. INPUT: Set input structure (.psf) and frame (.dcd, .pdb, .coor) files
# 3. INPUT: Set atom selections for RMSD, RMSF and RG
# 4. run with "./script.tcl"

## NOTE UNITS: RMSD, RMSF and RG values are in Angstrom (Å)


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


set COL_NAME_FRAME		"FRAME";	# Frame index
set COL_NAME_RESIDUE	"RESIDUE";	# Residue index
set COL_NAME_RMSD		"RMSD";		# RMSD (Å)
set COL_NAME_RMSF		"RMSF";		# RMSF (Å)
set COL_NAME_RG			"RG";		# RG (Å)


# =======================================
# INPUT
# =======================================
set psf_file		"../../common/amyld_wb.psf";		# TODO: input strcuture file (.psf)

# TODO: LIST of trajectory (.dcd) or single frame (.pdb, .coor) files separated by space	
#set frame_files		{ "../dna_wb_eq.dcd" "../dna_wb_eq1.dcd" };
set frame_files	[find_files ".." "amyld_wb_eq" ".dcd"];  # <dir> <prefix> <suffix> [min_num] [max_num]

set frame_index_ref		0;		# Reference Frame Index (default 0)
set frame_index_start 	-1;		# Inclusive, -1 for None
set frame_index_end 	-1;		# Exclusive, -1 for None

# TODO: ATOM SELECTION for RMSD, RG and RMSF
# Examples:
# 	1. protein and backbone and noh
# 	2. nucleic and backbone and noh
# 	3. protein and backbone and noh and not (resid 10 to 16)
set selection_rmsd		"protein and backbone and noh";		# TODO
set selection_rg		"$selection_rmsd";					# TODO

# TODO: RMSF BASE SELECTION for each Frame.
# -> Used to select the part of system conating residues
# -> Examples: "protein", "nucleic", "lipid", "all not water"
set rmsf_base_selection		"protein";						# TODO

# TODO: RMSF RESIDUE SELECTION within BASE_SELECTION
# -> Used to get RESIDUE_ID's of residues within BASE_SELECTION for which RMSF is required
# -> EXAMPLES
#    1. for Protein all residues: "protein and alpha" or "protein and name CA"
#    2. for Nucleic all residues: "nucleic and name C1'"
#    3. for RESIDUE RANGE, append: "and (resid <start> to <end>)" or "and (resid <one> <two> <three> ...)"
#       Example: "protein and alpha and (resid 10 to 14)", "nucleic and name C1' and not (resid 11 15)"
set rmsf_residue_selection 	"${rmsf_base_selection} and name CA"; # TODO


# ===========================================
# OUTPUT
# ===========================================
set out_file_prefix		"";

# Format
set out_format			"%.5f";		# RMSD Format.
set out_delimiter 		" ";		# output delimiter

# Plotting
set plot_output			1;			# [0/1] plot output using "plot_data.py"
set out_plot_file_type 	".pdf";		# [Optional] Only if plot_output = 1

set comment_token 		"#";		# Token used for Comments
set comment_header		0;			# Whether to comment out the columns header

set progress_report_frames  10000;	# report progress every N frames



# ==============================================
# MAIN
# ==============================================

# Output File names
if {[string trim $out_file_prefix] ne ""} {
	set out_file_prefix "[string trim $out_file_prefix].";
} else {
	set out_file_prefix "";
}

set out_file_rmsd "${out_file_prefix}rmsd.csv"
set out_file_rmsf "${out_file_prefix}rmsf.csv"
set out_file_rg "${out_file_prefix}rg.csv"


puts "\n=========================================================="
puts "=============  RMSD, RMSF and RG (bigdcd)  ================="
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

# All Atoms Selection -----------
set sel_all [atomselect $mol_id "all"];
set natom_all [$sel_all num];
if { $natom_all == 0 } {
	puts "\n-------------------------------------------------------------"
	puts " => ERROR: Molecule has no atoms : \"$psf_file\""
	puts "-------------------------------------------------------------\n"
	exit;
}

# RMSD and RG selections ----------------------------
set sel_rmsd [atomselect $mol_id "$selection_rmsd"];
set sel_rg [atomselect $mol_id "$selection_rg"];

set natom_rmsd [$sel_rmsd num];
set natom_rg [$sel_rg num];
if { $natom_rmsd == 0 } {
	puts "\n-------------------------------------------------------------"
	puts " => ERROR: RMSD Selection has no atoms : \"$selection_rmsd\""
	puts "-------------------------------------------------------------\n"
	exit;
}

if { $natom_rg == 0 } {
	puts "\n-------------------------------------------------------------"
	puts " => ERROR: RG Selection has no atoms : \"$selection_rg\""
	puts "-------------------------------------------------------------\n"
	exit;
}

# RMSF Selections ----------------------------
set sel_rmsf_base [atomselect $mol_id "$rmsf_base_selection"];
set natom_rmsf_base [$sel_rmsf_base num];
if { $natom_rmsf_base == 0 } {
	puts "\n-------------------------------------------------------------"
	puts " => ERROR: RMSF BASE_SELECTION has no atoms : \"$rmsf_base_selection\""
	puts "-------------------------------------------------------------\n"
	exit;
}

set res_ids [[atomselect $mol_id "$rmsf_residue_selection"] get resid];
set res_count [llength $res_ids];

if { $res_count == 0 } {
	puts "\n-------------------------------------------------------------"
	puts " => ERROR: RMSF RESIDUE_SELECTION has no residues : \"$rmsf_residue_selection\""
	puts "-------------------------------------------------------------\n"
	exit;
}





# =================================
# Output File
# =================================

# log to file
proc logf { f msg } {
	puts $f $msg;
}

# log to both console and file
proc log { f msg } {
	puts $msg; flush stdout;		# to stdout
	logf $f $msg;					# to output file
}

proc cleanup { } {
	global mol_id;
	global found_ref_frame ref_frame_filename;
	global file_rmsd file_rmsf file_rg;

	if { [info exists file_rmsd] == 1 } {
		flush $file_rmsd;
		close $file_rmsd;
	}

	if { [info exists file_rmsf] == 1 } {
		flush $file_rmsf;
		close $file_rmsf;
	}

	if { [info exists file_rg] == 1 } {
		flush $file_rg;
		close $file_rg;
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
set ref_frame_filename "ref_frame.temp-[clock seconds].coor";
set found_ref_frame 0;

# Bigdcd dryrun to find and save ref frame
proc try_save_ref_frame_bigdcd { i } {
	global frame_index_ref;

	set i [expr $i - 1];

	if { $i == $frame_index_ref } {
		global mol_id sel_all;
		global found_ref_frame ref_frame_filename;

		# Saving ref frame to pdb file
# 		set alla [atomselect $mol_id "all" frame $i];
		$sel_all writenamdbin $ref_frame_filename;
		set found_ref_frame 1;

		puts "\n-------------------------------------------------------------"
		puts " => bigdcd DRY-RUN: Reference Frame Found at Index : ${i}"
		puts "-------------------------------------------------------------\n"
		error "--";		# Exit bigdcd dryrun
	}
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

# =====================================
# PREPARE for MAIN-RUN
# =====================================
mol delete $mol_id;

set mol_id [mol new $psf_file];
mol addfile "$ref_frame_filename" type namdbin waitfor all molid $mol_id;
set sel_all [atomselect $mol_id "all"];

set sel_rg [atomselect $mol_id "$selection_rg"];

set sel_rmsd_ref [atomselect $mol_id "$selection_rmsd" frame 0];
set sel_rmsd [atomselect $mol_id "$selection_rmsd"];

set sel_rmsf_base_ref [atomselect $mol_id "$rmsf_base_selection" frame 0];
set sel_rmsf_base [atomselect $mol_id "$rmsf_base_selection"];

# Initialize RMSF Arrays
foreach r $res_ids {
	set rmsf_arr($r) 0;	# Output RMSF array

	set sel_res_str "${rmsf_base_selection} and resid $r and noh";
	set sel_rmsf_ref_res_arr($r) [atomselect $mol_id $sel_res_str frame 0];
	set sel_rmsf_res_arr($r) [atomselect $mol_id $sel_res_str];
}

set undo_rmsf_move 1;
if {[string match -nocase "$rmsf_base_selection*" "${selection_rmsd}"]} {
	puts "=> UNDO RMSF MOVE DISABLED"
    set undo_rmsf_move 0;
}

set frames_processed 0;

# -----------------------------------------
# Main function using bigdcd
# -----------------------------------------
proc calc_main { i } {
	global frame_index_ref frame_index_start frame_index_end;
	global mol_id sel_all;
	global sel_rmsd_ref sel_rmsd;
	global sel_rg;
	global sel_rmsf_base_ref sel_rmsf_base;
	global res_ids rmsf_arr sel_rmsf_ref_res_arr sel_rmsf_res_arr;
	global undo_rmsf_move;
	global file_rmsd file_rmsf file_rg;
	global out_delimiter out_format progress_report_frames;
	global frames_processed;

	set i [expr $i - 1];

	if { $i < $frame_index_start } {
		return;
	}

	if { $frame_index_end > 0 && $i >= $frame_index_end } {
		# Finished. Exit bigdcd
		error "Finished: FRAME RANGE \[$frame_index_start, $frame_index_end) processed";
	}

	# RG ------------------------------
	set rg_val [measure rgyr $sel_rg weight mass]
	set rg_str [format "$out_format" $rg_val]
	logf $file_rg "${i}${out_delimiter}${rg_str}";

	# RMSF ----------------------------
	set rmsf_transform [measure fit $sel_rmsf_base $sel_rmsf_base_ref];
	$sel_rmsf_base move $rmsf_transform;
	foreach r $res_ids {
		# compute the contribution to RMSF (loop through all residues)
		set sel_res_ref $sel_rmsf_ref_res_arr($r);
		set sel_res_cur $sel_rmsf_res_arr($r);

		set rmsd_r [measure rmsd $sel_res_cur $sel_res_ref]
		set rmsf_arr($r) [ expr $rmsf_arr($r) + ($rmsd_r * $rmsd_r) ]
	}

	# undo this move before calculating RMSD
	if { $undo_rmsf_move } {
		$sel_rmsf_base move [measure inverse $rmsf_transform];
	}
	unset rmsf_transform;

	# RMSD ----------------------------
	$sel_rmsd move [measure fit $sel_rmsd $sel_rmsd_ref];
	set rmsd_val [measure rmsd $sel_rmsd $sel_rmsd_ref];
	set rmsd_str [format "$out_format" $rmsd_val];
	logf $file_rmsd "${i}${out_delimiter}${rmsd_str}";


	# Post processing
	if {[expr $i % $progress_report_frames] == 0} {
		puts "\n---------------------------------------------------"
		puts " => processing Frame $i";
		puts "    RMSD: ${rmsd_str} Å | RG: ${rg_str} Å"
		puts "---------------------------------------------------\n"
	}

	incr frames_processed;
}


# ================================================
# OUTPUT File
# ================================================
set file_rmsd [open $out_file_rmsd w];
set file_rmsf [open $out_file_rmsf w];
set file_rg [open $out_file_rg w];

# ------------------------------------
# RMSD File header
# ------------------------------------
set comments_rmsd [list \
	"===================== RMSD (bigdcd) ======================="	\
    "INPUT Structure File  : \"${psf_file}\"" \
    "INPUT Frame File(s)   : \[${frame_files}\]" \
    "-----------------------------------------------------------" \
    "INPUT Atom Selection  : \"${selection_rmsd}\"" \
    "=> ATOM COUNT         : $natom_all (total), $natom_rmsd (selection)" \
    "-----------------------------------------------------------" \
	"INPUT REFERENCE FRAME Index : ${frame_index_ref}" \
	"INPUT FRAME RANGE           : \[${frame_index_start}, ${frame_end_str})" \
	"INPUT FRAME COUNT           : $frame_count_str" \
	"-----------------------------------------------------------" \
	"Units: RMSD (Angstrom Å)" \
	"-----------------------------------------------------------" \
]

foreach c $comments_rmsd {
	logf $file_rmsd "${comment_token} ${c}"
}

set header_rmsd "${COL_NAME_FRAME}${out_delimiter}${COL_NAME_RMSD}";
if {$comment_header == 1} {
	set header_rmsd "${comment_token}${header_rmsd}";
}

logf $file_rmsd $header_rmsd;


# ------------------------------------
# RG File header
# ------------------------------------
set comments_rg [list \
	"====================== RG (bigdcd) ========================"	\
    "INPUT Structure File  : \"${psf_file}\"" \
    "INPUT Frame File(s)   : \[${frame_files}\]" \
    "-----------------------------------------------------------" \
    "INPUT Atom Selection  : \"${selection_rg}\"" \
    "=> ATOM COUNT         : $natom_all (total), $natom_rg (selection)" \
    "-----------------------------------------------------------" \
	"INPUT FRAME RANGE     : \[${frame_index_start}, ${frame_end_str})" \
	"INPUT FRAME COUNT     : $frame_count_str" \
	"-----------------------------------------------------------" \
	"Units: RG (Angstrom Å)" \
	"-----------------------------------------------------------" \
]

foreach c $comments_rg {
	logf $file_rg "${comment_token} ${c}"
}

set header_rg "${COL_NAME_FRAME}${out_delimiter}${COL_NAME_RG}";
if {$comment_header == 1} {
	set header_rg "${comment_token}${header_rg}";
}

logf $file_rg $header_rg;




# ===========================================================
# MAIN-RUN: Calculate RMSD, RG, RMSF
# ===========================================================
eval "bigdcd calc_main auto [join $frame_files]";
bigdcd_wait;

if { $frames_processed == 0 } {
	puts "\n-------------------------------------------------------------"
	puts " => ERROR: No Frame Processed"
	puts "-------------------------------------------------------------\n"
	cleanup;
	exit;
}

# Final RMSF Calculation ------------
set avg_rmsf 0;

foreach r $res_ids {
    set rmsf_arr($r) [ expr sqrt($rmsf_arr($r) / $frames_processed) ]
    set avg_rmsf [expr $avg_rmsf + $rmsf_arr($r)]

    # Set the "user" field (of atoms of the residue) to RMSF value
#     set res_atoms $sel_rmsf_ref_res_arr($r)
#     $res_atoms set user $rmsf_arr($r);
}

set avg_rmsf [ expr $avg_rmsf / $res_count ]
set avg_rmsf_str "[format ${out_format} ${avg_rmsf}]"

# ------------------------------------
# RMSF File header
# ------------------------------------
set comments_rmsf [list \
	"===================== RMSF (bigdcd) ======================="	\
    "INPUT Structure File  : \"${psf_file}\"" \
    "INPUT Frame File(s)   : \[${frame_files}\]" \
    "-----------------------------------------------------------" \
    "INPUT Base Selection     : \"${rmsf_base_selection}\"" \
    "INPUT Residue Selection  : \"${rmsf_residue_selection}\"" \
    "=> ATOM Count            : $natom_all (total), $natom_rmsf_base (base selection)" \
    "=> RESIDUE Count         : $res_count" \
    "-----------------------------------------------------------" \
	"INPUT REFERENCE FRAME Index : ${frame_index_ref}" \
	"INPUT FRAME RANGE           : \[${frame_index_start}, ${frame_end_str})" \
	"INPUT FRAME COUNT           : $frame_count_str" \
	"=> FRAMES Processed         : ${frames_processed}" \
	"-----------------------------------------------------------" \
	"Units: RMSF (Angstrom Å)" \
	"=> Average RMSF: $avg_rmsf_str Å" \
	"-----------------------------------------------------------" \
]

foreach c $comments_rmsf {
	logf $file_rmsf "${comment_token} ${c}"
}

set header_rmsf "${COL_NAME_RESIDUE}${out_delimiter}${COL_NAME_RMSF}";
if {$comment_header == 1} {
	set header_rmsf "${comment_token}${header_rmsf}";
}

logf $file_rmsf $header_rmsf;

# RMSF DATA
foreach r $res_ids {
	logf $file_rmsf "${r}${out_delimiter}[format ${out_format} $rmsf_arr($r)]"
}



# GC ----------
cleanup;


# Plot Output -------------------

proc get_plot_filename {data_filename} {
	global out_plot_file_type
    return "[file rootname $data_filename]$out_plot_file_type"
}

if { $plot_output == 1 } {
	# RMSD Plot
	set plot_file_rmsd [get_plot_filename $out_file_rmsd];
	set cmd "./plot_data.py -t {RMSD Plot} -xl {Frame} -yl {RMSD (Å)} -ni -o \"${plot_file_rmsd}\" \"${out_file_rmsd}\"";
	eval $cmd;

	# RMSF Plot
	set plot_file_rmsf [get_plot_filename $out_file_rmsf];
	set cmd "./plot_data.py -t {RMSF Plot} -xl {Residue} -yl {RMSF (Å)} -ni -o \"${plot_file_rmsf}\" \"${out_file_rmsf}\"";
	eval $cmd;

	# RG Plot
	set plot_file_rg [get_plot_filename $out_file_rg];
	set cmd "./plot_data.py -t {RG Plot} -xl {Frame} -yl {RG (Å)} -ni -o \"${plot_file_rg}\" \"${out_file_rg}\"";
	eval $cmd;
}

set time_end [clock seconds];

puts "\n==================  FINISHED  ====================="
puts "OUTPUT DATA Files"
puts " ├─ RMSD : ${out_file_rmsd}"
puts " ├─ RMSF : ${out_file_rmsf}"
puts " └─ RG   : ${out_file_rg}"
if { $plot_output == 1 } {
	puts ""
	puts "OUTPUT PLOT Files"
	puts " ├─ RMSD : [get_plot_filename ${out_file_rmsd}]"
	puts " ├─ RMSF : [get_plot_filename ${out_file_rmsf}]"
	puts " └─ RG   : [get_plot_filename ${out_file_rg}]"
}

puts ""
puts "=> Time Taken: [expr $time_end -$time_start] secs"
puts "===================================================\n"

exit;
