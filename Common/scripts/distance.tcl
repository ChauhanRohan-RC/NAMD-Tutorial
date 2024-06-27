#############################################################################################
### Script to calculate DISTANCE B/W 2 ATOMS   from trajectory or coordinate file		   ##
#############################################################################################

## USAGE ------------------------
# 1. Copy script to working dir
# 2. INPUT: Set input strcuture (.psf) and frame (.dcd, .pdb, .coor) files
# 3. INPUT: Set selections for atom1 and atom2
# 4. run with "vmd -dispdev text -e distance.tcl"
# 5. OUTPUT: generates file "dist_vs_frame.dat"

## NOTE: Distance is in Angstrom (Å)

# ======================= INPUT ===========================
set psf_file		"dna_dummy.psf";		# input strcuture file (.psf)
# LIST of trajectory (.dcd) or single frame (.pdb, .coor) files separated by space
#set frame_files	{ "dna_gbis_pcf1.dcd" "dna_gbis_pcf2.dcd" };		
set frame_files	{ "dna_gbis_pcf.dcd" };

set selection_atom1	"resid 1 and name C5'";		# selection for first atom: SMD
set selection_atom2	"resid 16 and name C3'";	# selection for second atom: FIXED

### Direction Vector onto which ATOM 1 -> ATOM 2 linking vector is projected (dot prod)
## -> Must be a 3D VECTOR. 
## -> Comment out for no projection (abs distance between 1 and 2)

#set dir_vec 	{ 0.31759 0.88123 0.35008 };
set force_positive_dist		1;		# [Only when dir_vec is set] 1 : take only magntiude of distances. 0 : Distances can be negative

# ======================= OUTPUT ===========================
set out_file_name 		"dist_vs_frame.dat";	# output file name
set out_delimiter 		" \t ";					# output delimiter 

set comment_token 		"#";	# Token used for Comments
set comment_header		1;		# Whether to comment out the columns header
# -----------------------------------------------------

puts "==============  DISTANCE B/W ATOMS  ================="

# loading INPUT .psf and .dcd files
set mol_id [mol new $psf_file];	

foreach frame_file $frame_files {
	puts "---------------------------------------------"
	puts " -> Loading Frames from file : ${frame_file}"
	puts "---------------------------------------------"
	mol addfile $frame_file waitfor all molid $mol_id;
}


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

if {$nf == 1} {
	set start_frame_index 	0;
	set end_frame_index 	0;
} else {
	puts "----------------- Frame Indices ------------------";
	puts "INFO: Total Frames: $nf";
	puts "Frame index must be in range \[0, ${last_index}\], or -ve for back indices";
	puts ""
	puts -nonewline " -> START Frame Index (default: 0): "; flush stdout; set start_frame_index [gets stdin]
	puts -nonewline " -> END Frame Index (default: -1): "; flush stdout; set end_frame_index [gets stdin]
	puts ""
	
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
} 


# =========================  MAIN  ==================================

# open output file ofor writing
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

puts "--------------------------"
log2file "${comment_token}================  DISTANCE B/W 2 ATOMS  ================="
log "${comment_token}LOG: INPUT Structure File: \"${psf_file}\" | Frame File(s): \[${frame_files}\]"
log "${comment_token}LOG: ATOM1 SELECTION: \"${selection_atom1}\""
log "${comment_token}LOG: ATOM2 SELECTION: \"${selection_atom2}\""

if {[info exists dir_vec] && [llength $dir_vec] == 3 } {
	set dir_unit_vec [vecnorm $dir_vec]
	set has_dir_vec 1
	log "${comment_token}LOG: DIR Unit Vector: {$dir_unit_vec}  |  Force Positive Distances: ${force_positive_dist}"
} else {
	set dir_unit_vec { };	# Empty list
	set has_dir_vec 0
}

log "${comment_token}--------------------------"
log "${comment_token}LOG: Total Frames: ${nf}";
log "${comment_token}LOG: START Frame Index: ${start_frame_index} | END Frame Index: ${end_frame_index} | Frames for Calculation: [expr ${end_frame_index} - ${start_frame_index} + 1]"
log "${comment_token}--------------------------"

# Selecting Atoms
set atom1 [atomselect $mol_id $selection_atom1];
set atom2 [atomselect $mol_id $selection_atom2];

# Header for output file
if {$comment_header == 0} {
	set out_header "Frame${out_delimiter}DIST";
} else {
	set out_header "${comment_token}Frame${out_delimiter}DIST";
}

log2file $out_header;
#puts $out_header;

# Main Loop
for { set i $start_frame_index } { $i <= $end_frame_index  } { incr i } {
	$atom1 frame $i
	$atom2 frame $i
	
	set a1pos [lindex [$atom1 get {x y z}] 0]
	set a2pos [lindex [$atom2 get {x y z}] 0]
	
	set v [vecsub $a2pos $a1pos];		# ATOM 1 -> ATOM 2 Link Vector
	if {$has_dir_vec == 1} {
		set dist [vecdot $v $dir_unit_vec];		# Distance in Å
		if { $force_positive_dist == 1 } { 
			set dist [expr abs($dist)]
		}
	} else {
		set dsq [vecdot $v $v];				# Distance Square (self dot)
		set dist [expr { sqrt($dsq) }];		# Distance in Å
	}
	
	if {[expr $i % 10000] == 0} {
		puts "INFO: processing Frames $i-[expr min($i + 10000, $end_frame_index)]";
		
		if { $has_dir_vec == 1 } {
			puts "INFO: FRAME ${i} => 1 -> 2 Link Vector: {$v} | DIR Unit Vector: {$dir_unit_vec}"
			if {$force_positive_dist == 1} { set tok " | positive" } else { set tok "" }
			puts "INFO: FRAME ${i} => 1-2 DISTANCE (projected${tok}): ${dist} Å"
		} else {
			puts "INFO: FRAME ${i} => 1 -> 2 Link Unit Vector: {[vecnorm $v]}"
			puts "INFO: FRAME ${i} => 1-2 DISTANCE: ${dist} Å"
		}
	}
	
	set line "${i}${out_delimiter}${dist}"
	log2file $line
}

puts "=================  FINISHED  ===================="
puts "LOG: Output File: ${out_file_name}, delimiter: '$out_delimiter', comment token: '${comment_token}"
puts "================================================="

exit

