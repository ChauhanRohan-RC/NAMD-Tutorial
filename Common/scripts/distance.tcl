#!/usr/bin/env -S vmd -dispdev text -e

##############################################################################################
## Script to calculate DISTANCE B/W 2 ATOMS (Group of Atoms) from traj files				##
##############################################################################################
# **REQUIRE bigdcd.tcl**
package require bigdcd;
# OR
# source bigdcd.tcl;

# -> Calculates 1-2 interatomic Distance for each Frame
# -> Capture Checkpoints i.e. Frames with predefined distances of importance
# In case of group of atoms, Center of Mass / Geometric Center is used for distance calculation

## USAGE ------------------------
# 1. Copy script to working dir
# 2. INPUT: Set input strcuture (.psf) and frame (.dcd, .pdb, .coor) files
# 3. INPUT: Set selection 1 and 2
# 4. run with "./distance.tcl"
# 5. OUTPUT: generates file "dist_vs_frame.csv"
# 6. OUTPUT (optional) generates "checkpoint-<i>.pdb" file(s)

# Units: Distance => Angstrom (Å)


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


set COL_NAME_FRAME		"FRAME";		# Frame index
set COL_NAME_DIST		"DIST";		# Distance b/w atoms/group of atoms (Extension)



# =========================================
# INPUT
# =========================================
set psf_file		"../../common/amyld_wb.psf";		# TODO: input strcuture file (.psf)

# TODO: LIST of trajectory (.dcd) or single frame (.pdb, .coor) files separated by space
#set frame_files	 { "amyld_wb.pdb" "amyld_wb_eq.dcd" };
set frame_files	[find_files ".." "amyld_wb_eq" ".dcd"];  # <dir> <prefix> <suffix> [min_num] [max_num]

set frame_index_start 	-1;		# INCLUSIVE, -1 for None
set frame_index_end 	-1;		# Exclusive, -1 for None

set selection1	"protein and resid 1 and name CA";		# TODO: selection-1 [usually FIXED Atom(s)]
set selection2	"protein and resid 2 and name CA";		# TODO: selection-2  [usually SMD Atom(s)]

# In case of group of atoms, Use Center of Mass (COM) instead of geometric center
set use_cen_mass			1;

### Direction Vector onto which SEL-1 -> SEL-2 linking vector is projected (dot prod)
## -> Must be a 3D VECTOR. 
## -> Comment out for no projection (abs distance between 1 and 2)

#set dir_vec 	{ -1.0 0.0 0.0 };
set force_positive_dist		0;		# [Only when dir_vec is set] 1 : take only magntiude of distances. 0 : Distances can be negative



# =========================================
# OUTPUT
# =========================================
set out_file_prefix 	"dist_vs_frame";		# TODO: output file prefix

# Formatting
set out_format			"%.6f";		# Output distance format
set out_delimiter 		" ";		# output delimiter

### Checkpoint Frames Capture
# First frame(s) found to have the given Selection 1-2 distance(s) (Å) are saved as .pdb file
set checkpoints			{ }; 		# checkpoint 1-2 DIstances (in Å) separated by space
set checkpoint_tolerance		0.01;	# checkpoint Tolerance (in Å)

# Plotting
set plot_output				1;				# [0/1] plot output using "plot_data.py"
set show_interactive_plot 	1; 				# [0/1]

set comment_token 		"#";		# Token used for Comments
set comment_header		0;		# Whether to comment the columns header (for XmGrace). DO NOT comment header for my python scripts




# ==============================================
# MAIN
# ==============================================

# preprocess
set out_data_file 				"${out_file_prefix}.csv";			# output data file name
set out_plot_file 				"${out_file_prefix}.pdf";			# output figure file name
set checkpoint_file_prefix 		"${out_file_prefix}.checkpoint";	# checkpoint file prefix


puts "\n====================================================="
puts "==========  DISTANCE B/W ATOMS/GROUPS  ================"
puts "=====================================================\n"

set time_start [clock seconds];

## CHECKS ======================================

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

#foreach frame_file $frame_files {
#	puts "---------------------------------------------"
#	puts " -> Loading Frames from file : ${frame_file}"
#	puts "---------------------------------------------"
#	mol addfile $frame_file waitfor all molid $mol_id;
#}


# Selection -----------
set allatoms [atomselect $mol_id "all"];
set num_atoms_all [$allatoms num];
if { $num_atoms_all == 0 } {
	puts "\n-------------------------------------------------------------"
	puts " => ERROR: Molecule has no atoms : \"$psf_file\""
	puts "-------------------------------------------------------------\n"
	exit;
}

set sel1 [atomselect $mol_id $selection1];
set sel2 [atomselect $mol_id $selection2];

set num_atoms_sel1 [$sel1 num];
set num_atoms_sel2 [$sel2 num];

if { $num_atoms_sel1 == 0 } {
	puts "\n-------------------------------------------------------------"
	puts " => ERROR: Selection 1 has no atoms : \"$selection1\""
	puts "-------------------------------------------------------------\n"
	exit;
}

if { $num_atoms_sel2 == 0 } {
	puts "\n-------------------------------------------------------------"
	puts " => ERROR: Selection 2 has no atoms : \"$selection2\""
	puts "-------------------------------------------------------------\n"
	exit;
}


# Direction Vector ---------
if {[info exists dir_vec] && [llength $dir_vec] == 3 } {
	set dir_unit_vec [vecnorm $dir_vec];
	set has_dir_vec 1;
} else {
	set dir_unit_vec { };	# Empty list
	set has_dir_vec 0;
}

# Checkpoints --------------
if {[info exists checkpoints] == 1 && [llength $checkpoints] > 0} {
	set checkpoints_copy  $checkpoints;
	set has_checkpoints 1;
} else {
	set checkpoints_copy [list];
	set has_checkpoints 0;
}


# OUTPUT Files ================================
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

puts "\n---------------------------------------------------------"
log "${comment_token} ============  DISTANCE B/W 2 ATOMS (or Group of Atoms)  ============="
log "${comment_token} INPUT Structure File: \"${psf_file}\""
log "${comment_token} INPUT Frame File(s): \[${frame_files}\]"
log "${comment_token} INPUT SELECTION-1 : \"${selection1}\""
log "${comment_token} INPUT SELECTION-2 : \"${selection2}\""
log "${comment_token} => ATOM COUNT: Total ($num_atoms_all), Selection-1 ($num_atoms_sel1), Selection-2 ($num_atoms_sel2)"
log "${comment_token}-------------------------------------------------------"
log "${comment_token} FRAME Index RANGE: \[${frame_index_start}, ${frame_end_str})  |  COUNT: $frame_count_str"
if { $has_dir_vec == 1 } {
	log "${comment_token} DIR Unit Vector: {$dir_unit_vec}  |  Force Positive Distances: ${force_positive_dist}"
}
log "${comment_token}-------------------------------------------------------"
log "${comment_token} OUTPUT Units: Distance (Angstrom Å)"
log "${comment_token}-------------------------------------------------------"
puts ""

# Header for output file
set out_header "${COL_NAME_FRAME}${out_delimiter}${COL_NAME_DIST}";
if {$comment_header == 1} {
	set out_header "${comment_token}${out_header}";	
} 

log2file $out_header;
#puts $out_header;


# =======================================
# Function to calculate Distance
# =======================================

proc calc_dist { i } {
	global frame_index_start frame_index_end;

	set i [expr $i -1];

	if { $i < $frame_index_start } {
		return;
	}

	if { $frame_index_end > 0 && $i >= $frame_index_end } {
		# Finished. Exit bigdcd
		error "Finished: FRAME RANGE \[$frame_index_start, $frame_index_end) processed";
	}


	global sel1 sel2 allatoms num_atoms_sel1 num_atoms_sel2 use_cen_mass;
	global has_dir_vec dir_unit_vec force_positive_dist;
	global has_checkpoints checkpoints checkpoints_copy checkpoint_tolerance checkpoint_file_prefix;
	global out_delimiter out_format;

	$sel1 frame $i
	$sel2 frame $i

	if { $num_atoms_sel1 > 1 && $use_cen_mass == 1 } {
		set sel1_pos [measure center $sel1 weight mass];	# COM
	} else {
		set sel1_pos [measure center $sel1];				# Geo Center
	}

	if { $num_atoms_sel2 > 1 && $use_cen_mass == 1 } {
		set sel2_pos [measure center $sel2 weight mass];	# COM
	} else {
		set sel2_pos [measure center $sel2];				# Geo Center
	}


	set v [vecsub $sel2_pos $sel1_pos];			# SEL-1 -> SEL-2 Link Vector
	if {$has_dir_vec == 1} {
		set dist [vecdot $v $dir_unit_vec];		# Distance in Å
		if { $force_positive_dist == 1 } {
			set dist [expr abs($dist)]
		}
	} else {
		set dsq [vecdot $v $v];				# Distance Square (self dot)
		set dist [expr { sqrt($dsq) }];		# Distance in Å
	}

	set dist_str "[format ${out_format} ${dist}]";


	if {$has_checkpoints == 1} {
		for { set ci 0 } { $ci < [llength $checkpoints_copy] } { } {
			set cp [lindex $checkpoints_copy $ci]
			if { [expr abs($dist - $cp)] <= $checkpoint_tolerance } {
				# we got the checkpoint, Save the frame
				set check_id [expr [lsearch $checkpoints $cp] + 1]

				#$allatoms frame $i;		# causes errors, not needed with bigdcd
				set check_file "${checkpoint_file_prefix}-${check_id}.pdb";
				$allatoms writepdb $check_file;
				puts "\n--------------------------------------------------------------"
				puts " -> CHECKPOINT ${check_id} Found => FRAME Index: $i | Output File: \"${check_file}\" | 1-2 DISTANCE: ${cp} Å (Requested), ${dist} Å (Actual)"
				puts "--------------------------------------------------------------\n"

				set checkpoints_copy [lreplace $checkpoints_copy $ci $ci];	# remove this checkpoint
			} else {
				incr ci;	# Increment
			}
		}
	}

	if {[expr $i % 10000] == 0} {
		puts "\n--------------------------------------------------------------"
		if { $has_dir_vec == 1 } {
			puts "=> FRAME ${i} => 1 -> 2 Link Vector: {$v} | DIR Unit Vector: {$dir_unit_vec}"
			if {$force_positive_dist == 1} { set tok " | positive" } else { set tok "" }
			puts "=> FRAME ${i} => 1-2 DISTANCE (projected${tok}): ${dist_str} Å"
		} else {
			puts "=> FRAME ${i} => 1 -> 2 Link Unit Vector: {[vecnorm $v]}"
			puts "=> FRAME ${i} => 1-2 DISTANCE: ${dist_str} Å"
		}
		puts "--------------------------------------------------------------\n"
	}

	set line "${i}${out_delimiter}${dist_str}";
	log2file $line;
}


# =====================================
# Execute Calculation
# =====================================
eval "bigdcd calc_dist auto [join $frame_files]";
bigdcd_wait;

# GC
flush $out_file;
close $out_file;
mol delete $mol_id;


# Plot Output -------------------
if { $plot_output == 1 } {
	set cmd "./plot_data.py -x ${COL_NAME_FRAME} -y ${COL_NAME_DIST}  -xl {Frame} -yl {Dist (Å)} -o \"${out_plot_file}\" \"${out_data_file}\"";
	if { $show_interactive_plot == 0 } {
		set cmd "$cmd -ni";
	}

	eval $cmd;
}


set time_end [clock seconds];

puts "\n=======================  FINISHED  =========================="
puts "=> Output File: ${out_data_file}"
if { $plot_output == 1 } {
	puts "=> OUTPUT Plot File: \"${out_plot_file}\""
}
puts "-> Time Taken: [expr $time_end -$time_start] secs"
puts "=============================================================\n"

exit;

