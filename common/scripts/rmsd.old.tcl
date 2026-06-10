#!/usr/bin/env -S vmd -dispdev text -e

#############################################################################
# Script to calculate RMSD of a molecule from trajectory in VMD			   ##
#############################################################################
# NOTE: Plotting data requires "plot_data.py"

### USAGE ------------------
# 1. Copy this script in working dir
# 2. INPUT: Set input structure (.psf) and frame (.dcd, .pdb, .coor) files
# 3. INPUT: Set ATOM_SELECTION for each frame
# 4. run with "./rmsd.tcl"
# 5. Generates output file "rmsd.dat"

## NOTE UNITS: RMSD values are in Angstrom (Å)
proc find_files { directory prefix suffix { sort_natural 1 } {return_abs_path 0} } {
    if { $return_abs_path == 1 } { set directory [file normalize $directory]; }
	set file_list [glob -nocomplain -join "${directory}" "${prefix}*${suffix}"];
    if { $sort_natural == 1 } { set file_list [lsort -dictionary $file_list] };
    return $file_list;
}


# ======================= INPUT ===========================
set psf_file		"amyld_wb.psf";		# TODO: input strcuture file (.psf)

# TODO: LIST of trajectory (.dcd) or single frame (.pdb, .coor) files separated by space	
#set frame_files		{ "../dna_wb_eq.dcd" "../dna_wb_eq1.dcd" };
set frame_files	[find_files "." "amyld_wb_eq" ".dcd"];  # <dir> <prefix> <suffix>

# TODO: ATOM SELECTION for each Frame
set atom_selection		"protein and backbone and noh";
# Examples:
# 	1. protein and backbone and noh
# 	2. nucleic and backbone and noh
# 	3. protein and backbone and noh and not (resid 10 to 16)

# ======================= OUTPUT ===========================
set out_file_prefix		"rmsd"
set out_data_file 		"${out_file_prefix}.csv";		# output data file

# Plotting
set plot_output				1;				# [0/1] plot output using "plot_data.py"
set show_interactive_plot 	0; 				# [0/1]
set out_plot_file 			"${out_file_prefix}.pdf";	# [Optional] Only if plot_output = 1

set out_delimiter 		" ";			# output delimiter
set comment_token 		"#";	# Token used for Comments
set comment_header		0;		# Whether to comment out the columns header
# -----------------------------------------------------

puts "=============  RMSD  ================="

#puts -nonewline " -> set Molecule ID (default: top): "; flush stdout; 
#set mol_id [gets stdin]
#if { [string trim $mol_id] eq ""} { set mol_id top }
#puts ""

# loading INPUT .psf and .dcd files
set mol_id [mol new $psf_file];	

foreach frame_file $frame_files {
	puts "---------------------------------------------"
	puts " -> Loading Frames from file : ${frame_file}"
	puts "---------------------------------------------"
	mol addfile $frame_file waitfor all molid $mol_id;
}

#puts "---------------- ATOM SELECTION ---------------"
#puts " Examples: "
#puts " \t 1. protein and backbone and noh"
#puts " \t 2. nucleic and backbone and noh"
#puts " \t 3. protein and backbone and noh and not (resid 10 to 76)"
#puts ""
#puts -nonewline " -> set ATOM SELECTION for each frame: "; flush stdout; 
#set atom_selection [gets stdin]
#puts ""

### FRAMES INPUT ---------------------------
# Number of frames
set nf [molinfo $mol_id get numframes]
set last_index [expr $nf - 1]

if {$nf == 0} {
	puts "-------------------------------------------"
	puts "INFO: There are no frames in input file(s): \[${frame_files}\]"
	puts "LOG: QUITING (NO FRAMES)"
	puts "-------------------------------------------"
	exit
}

puts "----------------- Frame Indices ------------------";
puts "INFO: Total Frames: $nf";
puts "Frame index must be in range \[0, ${last_index}\], or -ve for back indices";
puts ""
puts -nonewline " -> Reference Frame Index (default: 0): "; flush stdout; set ref_frame_index [gets stdin]
puts -nonewline " -> START Frame Index (default: 0): "; flush stdout; set start_frame_index [gets stdin]
puts -nonewline " -> END Frame Index (default: -1): "; flush stdout; set end_frame_index [gets stdin]
puts ""

# NORMALIZING REFERENCE_FRAME_INDEX (default 0)
if { [string trim $ref_frame_index] eq ""} {
	set ref_frame_index 0
} elseif {$ref_frame_index < 0} {
	set ref_frame_index [expr $last_index + (($ref_frame_index + 1) % -$nf)]
} elseif { $ref_frame_index >= $nf } {
	set ref_frame_index [expr $ref_frame_index % $nf]
}

# NORMALIZING START_FRAME_INDEX (default 0)
if { [string trim $start_frame_index] eq ""} {
	set start_frame_index 0
} elseif {$start_frame_index < 0} {
	set start_frame_index [expr $last_index + (($start_frame_index + 1) % -$nf)]
} elseif { $start_frame_index >= $nf } {
	set start_frame_index [expr $start_frame_index % $nf]
}

# Normalizing END_FRAME_INDEX (default last_index)
if { [string trim $end_frame_index] eq ""} {
	set end_frame_index $last_index
} elseif {$end_frame_index < 0} {
	set end_frame_index [expr $last_index + (($end_frame_index + 1) % -$nf)]
} elseif { $end_frame_index >= $nf } {
	set end_frame_index [expr $end_frame_index % $nf]
}

# =========================  MAIN  ==================================

# open output file ofor writing
set out_file [open $out_data_file w]

proc log2file { msg } { 
	global out_file;
	puts $out_file $msg;	# to output file
}

proc log { msg } { 
	global out_file;
	puts $msg; flush stdout;		# to stdout
	log2file $msg;					# to output file
}

log2file "${comment_token}================ RMSD =================="
puts "--------------------------"
log "${comment_token}LOG: INPUT Structure File: \"${psf_file}\" | Frame File(s): \[${frame_files}\]"
log "${comment_token}LOG: ATOM SELECTION: \"${atom_selection}\""
log "${comment_token}--------------------------"
log "${comment_token}LOG: Total Frames: ${nf}";
log "${comment_token}LOG: REFERENCE Frame Index: ${ref_frame_index}";
log "${comment_token}LOG: START Frame Index: ${start_frame_index} | END Frame Index: ${end_frame_index} | Frames for RMSD: [expr ${end_frame_index} - ${start_frame_index} + 1]"

# Reference Conformation
set allatoms [atomselect $mol_id "all"];
set ref_frame [atomselect $mol_id $atom_selection frame $ref_frame_index];
log "${comment_token}LOG: selected ATOM COUNT in each frame: [$ref_frame num]"
log "${comment_token}--------------------------"

# current frame variable, must have same atom selection as reference
set cur_frame [atomselect $mol_id $atom_selection];

# Header for output file
if {$comment_header == 0} {
	set out_header "FRAME${out_delimiter}RMSD";
} else {
	set out_header "${comment_token}FRAME${out_delimiter}RMSD";
}

log2file $out_header;	
#puts $out_header;

# Main Loop
for { set i $start_frame_index } { $i <= $end_frame_index } { incr i } {
	$allatoms frame $i;
	$cur_frame frame $i;
	$allatoms move [measure fit $cur_frame $ref_frame];
	
	set rmsd_val [measure rmsd $cur_frame $ref_frame]
	
	if {[expr $i % 10000] == 0} {
		puts "INFO: processing Frames $i-[expr min($i + 10000, $end_frame_index)]";
		puts "INFO: FRAME ${i} -> RMSD: ${rmsd_val} Å"
	}
	
	set out_line "${i}${out_delimiter}${rmsd_val}";
	log2file $out_line;
	# puts out_line;
}

# GC
flush $out_file;
close $out_file;
mol delete $mol_id;

# Plot Output -------------------
if { $plot_output == 1 } {
	set cmd "./plot_data.py -t {RMSD Plot} -xl {Frame} -yl {RMSD (Å)} -o \"${out_plot_file}\" \"${out_data_file}\"";
	if { $show_interactive_plot == 0 } {
		set cmd "$cmd -ni";
	}

	eval $cmd;
}

puts "\n================  FINISHED  ==================="
puts "=> OUTPUT Data File: \"${out_data_file}\""
if { $plot_output == 1 } {
	puts "=> OUTPUT Plot File: \"${out_plot_file}\""
}
puts "===============================================\n"

exit;
