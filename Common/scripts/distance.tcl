#############################################################################################
### Script to calculate DISTANCE B/W 2 ATOMS   from trajectory or coordinate file		   ##
#############################################################################################

## USAGE ------------------------
# 1. Copy script to working dir
# 2. INPUT: Set input .psf and .dcd (or .pdb .coor coordinate) files
# 3. INPUT: Set selections for atom1 and atom2
# 4. run with "vmd -dispdev text -e distance.tcl"
# 5. OUTPUT: generates file "dist_vs_frame.dat"

# ======================= INPUT ===========================
set in_psf_file "../common/dna.psf";		# input strcuture file (.psf)
set in_frame_file "dna_gbis_pcf.dcd";		# input trajectory (.dcd) or a single frame (.pdb, .coor)

set selection_atom1 "nucleic and resid 1 and name C5'";		# selection for first atom
set selection_atom2 "nucleic and resid 16 and name C3'";	# selection for second atom

# ======================= OUTPUT ===========================
set out_file_name "dist_vs_frame.dat";		# output file name
set out_delimiter " \t ";					# output delimiter 
# -----------------------------------------------------



puts "#=============  DISTANCE B/W ATOMS  ================="

# loading INPUT .psf and .dcd files
set mol_id [mol new $in_psf_file];	

puts "---------------------------------------------"
puts " -> Loading Frames from file : ${in_frame_file}"
puts "---------------------------------------------"
mol addfile $in_frame_file waitfor all molid $mol_id;


### FRAMES INPUT ---------------------------
# Number of frames
set nf [molinfo $mol_id get numframes]
set last_index [expr $nf - 1]

if {$nf == 0} {
	puts "-------------------------------------------"
	puts "INFO: There are no frames in file: ${in_frame_file}"
	puts "LOG: QUITING (No Frames)"
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
log2file "#================  DISTANCE B/W 2 ATOMS  ================="
log "#LOG: INPUT Structure File: ${in_psf_file} | Frame(s) File: ${in_frame_file}$"
log "#--------------------------"
log "#LOG: Mol ID: $mol_id";
log "#LOG: ATOM1 SELECTION: \"${selection_atom1}\""
log "#LOG: ATOM2 SELECTION: \"${selection_atom2}\""
log "#--------------------------"
log "#LOG: Frame Count: ${nf}";
log "#LOG: START Frame Index: ${start_frame_index} | END Frame Index: ${end_frame_index} | Frames for Calculation: [expr ${end_frame_index} - ${start_frame_index} + 1]"
log "#--------------------------"

# Selecting Atoms
set atom1 [atomselect $mol_id $selection_atom1];
set atom2 [atomselect $mol_id $selection_atom2];

# Header for output file
set out_header "#Frame${out_delimiter}DIST";
log2file $out_header;
#puts $out_header;

# Main Loop
for { set i $start_frame_index } { $i <= $end_frame_index  } { incr i } {
	if {[expr $i % 10000] == 0} {
		puts "INFO: processing Frames $i-[expr min($i + 10000, $end_frame_index)]";
	}

	$atom1 frame $i
	$atom2 frame $i
	
	set a1pos [lindex [$atom1 get {x y z}] 0]
	set a2pos [lindex [$atom2 get {x y z}] 0]
	
	set v [vecsub $a2pos $a1pos]
	set dsq [vecdot $v $v]
	set dist [expr { sqrt($dsq) }]
	
	set line "${i}${out_delimiter}${dist}"
	log2file $line
}

puts "====================  FINISHED  ==========================="
exit

