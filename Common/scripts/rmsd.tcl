## Script to calculate RMSD of a molecule at each timestep (TS) in VMD

### USAGE ------------------
# 1. Copy this script in working dir
# 2. Load the .psf and .dcd in VMD
# 3. Go to Extensions > Tk Console
# 4. enter command: "source rmsd.tcl"
# 5. enter ATOM_SELECTION and FRAME indices
# 6. Generates output file "rmsd.dat"

# ======================= CONFIG ===========================
set out_file_name "rmsd.dat";
set out_delimiter " \t ";

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

log "#=============  RMSD  ================="

puts -nonewline " -> set Molecule ID (default: top): "; flush stdout; 
set mol_id [gets stdin]
if { [string trim $mol_id] eq ""} { set mol_id top }
puts ""

puts "---------------- ATOM SELECTION ---------------"
puts " Examples: "
puts " \t 1. protein and backbone and noh"
puts " \t 2. nucleic and backbone and noh"
puts " \t 3. protein and backbone and noh and not (resid 10 to 76)"
puts ""
puts -nonewline " -> set ATOM SELECTION for each frame: "; flush stdout; 
set atom_selection [gets stdin]
puts ""

### FRAMES INPUT ---------------------------
# Number of frames
set nf [molinfo $mol_id get numframes]
set last_index [expr $nf - 1]

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
puts "--------------------------"
log "#LOG: Mol ID: $mol_id";
log "#LOG: ATOM SELECTION: '${atom_selection}'"
log2file "#--------------------------"
log "#LOG: Frame Count: ${nf}";
log "#LOG: REFERENCE Frame Index: ${ref_frame_index}";
log "#LOG: START Frame Index: ${start_frame_index} | END Frame Index: ${end_frame_index} | Frames for RMSD: [expr ${end_frame_index} - ${start_frame_index} + 1]"
log2file "#--------------------------"

# Reference Conformation
set ref_frame [atomselect $mol_id $atom_selection frame $ref_frame_index];
log "#LOG: selected ATOM COUNT in each frame: [$ref_frame num]"

# current frame variable, must have same atom selection as reference
set cur_frame [atomselect $mol_id $atom_selection];

# Header for output file
set out_header "#Frame${out_delimiter}RMSD";
log2file $out_header;	
#puts $out_header;

# Main Loop
for { set i $start_frame_index } { $i <= $end_frame_index } { incr i } {
	$cur_frame frame $i;
	$cur_frame move [measure fit $cur_frame $ref_frame];
	
	set out_line "${i}${out_delimiter}[measure rmsd $cur_frame $ref_frame]";
	log2file $out_line;
	# puts out_line;
}

close $out_file;

puts "---------------------------------"
puts "LOG: FINISHED. Output File: ${out_file_name}, delimiter: '$out_delimiter'"
