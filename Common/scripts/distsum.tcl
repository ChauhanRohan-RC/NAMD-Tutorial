##################################################################################################
## Script to calculate SUM OF DISTANCES B/W Pair of ATOMS (Group of Atoms) from traj files		##
##################################################################################################
# **REQUIRE bigdcd.tcl**
package require bigdcd;
# OR
# source bigdcd.tcl;

# -> Calculates sum of interatomic (inter group) distances for each Frame
# -> Capture Checkpoints i.e. Frames with predefined distance_sum of importance

## WORKING -------------------
# 1. Corresponding elements in selection-1 and selection-2 lists form a selection pair
# 2. For each selection pair:
#	 -> Inter-atomic (or inter-group) distance is calculated
#		(for group of atoms, Center of Mass / Geometric Center is used for distance calculation)
#    -> [OPTIONAL] distance is projected on to dir_vec (if specified)
# 3. distance for all selection pairs, and the sum are reported for each frame

## USAGE ------------------------
# 1. Copy script to working dir
# 2. INPUT: Set input strcuture (.psf) and frame (.dcd, .pdb, .coor) files
# 3. INPUT: Set selection 1 and 2 lists
# 4. run with "vmd -dispdev text -e distance.tcl"
# -----------------
# 5. OUTPUT: generates file "distsum_vs_frame.csv"
# 6. OUTPUT (optional) generates "distsum_checkpoint-<i>.pdb" file(s)

# Units: Distance => Angstrom (Å)


set COL_NAME_FRAME			"FRAME";		# Frame index
set COL_NAME_DIST_PREFIX	"DIST";			# Distance b/w atoms/group of atoms (Extension)
set COL_NAME_DIST_SUM		"DIST_SUM";		# Distance Sum

# ======================= INPUT ===========================
set psf_file		"amyld_wb.psf";		# TODO: input strcuture file (.psf)

# TODO: LIST of trajectory (.dcd) or single frame (.pdb, .coor) files separated by space
set frame_files	{ "amyld_wb.pdb" "amyld_wb_eq.restart.coor" "amyld_wb_eq.dcd" };

set frame_index_start 	-1;		# INCLUSIVE, -1 for None
set frame_index_end 	-1;		# Exclusive, -1 for None

# TODO: Selection lists
# -> corresponding elements in sel-1 sel-2 lists form a pair, for which distances are calculated
set selection1_list	{
	"protein and (resid 1 and name N)"
	"protein and (resid 47 and name NH1)"
	"protein and (resid 70 and name N)"
};

set selection2_list	{
	"protein and (resid 65 and name OD2)"
	"protein and (resid 7 and name OD1)"
	"protein and (resid 25 and name O)"
};


# In case of group of atoms, Use Center of Mass (COM) instead of geometric center
set use_cen_mass			1;		# [0/1]

### Direction Vector onto which SEL-1(i) -> SEL-2(i) linking vector is projected (dot prod)
## -> Must be a 3D VECTOR.
## -> Comment out for no projection (abs distance between selection pairs)

#set dir_vec 	{ -1.0 0.0 0.0 };
set force_positive_dist		0;		# # [0/1] [Only with dir_vec] 1 : take only magntiude of distances for each selection pair. 0 : Distances can be negative

# ======================= OUTPUT ==============================================
set out_file_name 		"distsum_vs_frame.csv";		# TODO: output file name
set out_only_dist_sum	0;	 # [0/1] Output only the sum of distance, else also output individual distances

set out_delimiter 		" ";						# output delimiter

### Checkpoint Frames Capture
# First frame(s) found to have the given distance sum (Å) are saved as .pdb file
set checkpoints			{ }; 			# checkpoint dIstance sum (in Å) separated by space
set checkpoint_tolerance		0.01;	# checkpoint Tolerance (in Å)
set checkpoint_out_file_name_prefix 	"distsum_checkpoint";

set comment_token 		"#";	# Token used for Comments
set comment_header		0;		# [0/1] Whether to comment the columns header (for XmGrace). DO NOT comment header for my python scripts
# -----------------------------------------------------




# ==============================================
# MAIN
# ==============================================

puts "\n====================================================="
puts "==========  SUM of DISTANCE B/W ATOMS/GROUPS  ================"
puts "=====================================================\n"

## CHECKS ======================================

# Selection List
if { [llength $selection1_list] != [llength $selection2_list] || [llength $selection1_list] == 0 } {
	puts "\n------------------------------------------"
	puts " => Selection lists must be of equal length and non-empty. Given"
	puts "    -> Selection-1 List : \[${selection1_list}\]"
	puts "    -> Selection-2 List : \[${selection2_list}\]"
	puts "------------------------------------------\n"
	exit;
}

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
if { [info exists frame_index_end] == 0 || $frame_index_end <= 0 } {
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

set sel1l [list];		# List of sel-1
set sel2l [list];		# List of sel-2
foreach s1 $selection1_list s2 $selection2_list {
	set sel1 [atomselect $mol_id $s1];
	set sel2 [atomselect $mol_id $s2];

	if { [$sel1 num] == 0 } {
		puts "\n-------------------------------------------------------------"
		puts " => ERROR: Selection has no atoms : \"$s1\""
		puts "-------------------------------------------------------------\n"
		exit;
	}

	if { [$sel2 num] == 0 } {
		puts "\n-------------------------------------------------------------"
		puts " => ERROR: Selection has no atoms : \"$s2\""
		puts "-------------------------------------------------------------\n"
		exit;
	}

	lappend sel1l $sel1;
	lappend sel2l $sel2;

	unset sel1; unset sel2;
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
set out_file [open $out_file_name w]

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
log2file "${comment_token}============  DISTANCE SUM  ============="
log "${comment_token} INPUT Structure File: \"${psf_file}\""
log "${comment_token} INPUT Frame File(s): \[${frame_files}\]"
log "${comment_token}-------------------------------------------------------"
log "${comment_token} NOTE: Corresponding elements in SELECTION-1 and SELECTION-2 lists form a selection pair"
log "${comment_token} INPUT SELECTION-1 List : \[${selection1_list}\]"
log "${comment_token} INPUT SELECTION-2 List : \[${selection2_list}\]"
log "${comment_token}-------------------------------------------------------"
log "${comment_token} FRAME Index RANGE: \[${frame_index_start}, ${frame_end_str})  |  COUNT: $frame_count_str"
if { $has_dir_vec == 1 } {
	log "${comment_token} DIR Unit Vector: {$dir_unit_vec}  |  Force Positive Distances: ${force_positive_dist}"
}
log "${comment_token} UNITS: Distance = Angstrom (Å)"
log "${comment_token}-------------------------------------------------------"
puts ""

# Header for output file
set lhead [list $COL_NAME_FRAME];
if { $out_only_dist_sum != 1 } {
	for {set i 1} { $i <= [llength $selection1_list] } {incr i} {
		lappend lhead "${COL_NAME_DIST_PREFIX}${i}";
	}
}

lappend lhead $COL_NAME_DIST_SUM;
set out_header [join $lhead "$out_delimiter"]

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

	if { $i < $frame_index_start || ($frame_index_end > 0 && $i >= $frame_index_end) } {
		return;
	}


	global sel1l sel2l allatoms use_cen_mass out_only_dist_sum;
	global has_dir_vec dir_unit_vec force_positive_dist;
	global has_checkpoints checkpoints checkpoints_copy checkpoint_tolerance checkpoint_out_file_name_prefix;
	global out_delimiter;

	set distl [list];
	set dist_sum 0;

	foreach sel1 $sel1l sel2 $sel2l {
		$sel1 frame $i;
		$sel2 frame $i;

		if { [$sel1 num] > 1 && $use_cen_mass == 1 } {
			set sel1_pos [measure center $sel1 weight mass];	# COM
		} else {
			set sel1_pos [measure center $sel1];				# Geo Center
		}

		if { [$sel2 num] > 1 && $use_cen_mass == 1 } {
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

		lappend distl $dist;
		set dist_sum [expr $dist_sum + $dist]
	}



	if {$has_checkpoints == 1} {
		for { set ci 0 } { $ci < [llength $checkpoints_copy] } { } {
			set cp [lindex $checkpoints_copy $ci]
			if { [expr abs($dist_sum - $cp)] <= $checkpoint_tolerance } {
				# we got the checkpoint, Save the frame
				set check_id [expr [lsearch $checkpoints $cp] + 1]

				#$allatoms frame $i;		# causes errors, not needed with bigdcd
				set check_file "${checkpoint_out_file_name_prefix}-${check_id}.pdb";
				$allatoms writepdb $check_file;
				puts "\n--------------------------------------------------------------"
				puts " -> CHECKPOINT ${check_id} Found => FRAME Index: $i | Output File: \"${check_file}\" | 1-2 DISTANCE: ${cp} Å (Requested), ${dist_sum} Å (Actual)"
				puts "--------------------------------------------------------------\n"

				set checkpoints_copy [lreplace $checkpoints_copy $ci $ci];	# remove this checkpoint
			} else {
				incr ci;	# Increment
			}
		}
	}

	if {[expr $i % 10000] == 0} {
		puts "INFO: processing Frames $i";

		if { $has_dir_vec == 1 } {
			if { $force_positive_dist == 1 } { set tok " | positive" } else { set tok "" }
			puts "INFO: FRAME ${i} => DISTANCE SUM (projected${tok}): ${dist_sum} Å"
		} else {
			puts "INFO: FRAME ${i} => DISTANCE SUM: ${dist_sum} Å"
		}
	}

	if { $out_only_dist_sum == 1 } {
		set line "${i}${out_delimiter}${dist_sum}";
	} else {
		set distl_str [join $distl "$out_delimiter"]

		set line "${i}${out_delimiter}${distl_str}${out_delimiter}${dist_sum}";
	}

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

puts "=================  FINISHED  ===================="
puts "LOG: Output File: ${out_file_name}, delimiter: '$out_delimiter', comment token: '${comment_token}"
puts "================================================="

exit;

