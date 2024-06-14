## Script to calculate RMSF of each residue of a molecule over each timestep (TS) in VMD
## NOTE: the function also sets the "user" field of each atom to the RMSF value for its residue (Like a TAG)

### USAGE-------------------
# 1. Copy this script in working dir
# 2. INPUT: Set input structure (.psf) and frame (.dcd, .pdb, .coor) files
# 3. INPUT: Set BASE_SELECTION and RESIDUE_SELECTION
# 4. run with "vmd -dispdev text -e rmsf.tcl"
# 5. Generates output file "rmsf.dat"

## NOTE: RMSF values are in Angstrom (Å)

# ======================= INPUT ===========================
set psf_file		"../common/dna_wb.psf";		# TODO: input strcuture file (.psf)

# TODO: LIST of trajectory (.dcd) or single frame (.pdb, .coor) files separated by space	
set frame_files		{ "dna_wb_eq.dcd" "dna_wb_eq1.dcd" "dna_wb_eq2.dcd" };	

# TODO: BASE SELECTION for each Frame.
# -> Used to select the part of system conating residues
# -> Examples: "protein", "nucleic", "lipid", "all not water"
set base_selection		"nucleic";

# TODO: RESIDUE SELECTION within BASE_SELECTION
# -> Used to get RESIDUE_ID's of residues within BASE_SELECTION for which RMSF is required
# -> EXAMPLES
#    1. for Protein all residues: "protein and alpha" or "protein and name CA" 
#    2. for Nucleic all residues: "nucleic and name C1'"
#    3. for RESIDUE RANGE, append: "and (resid <start> to <end>)" or "and (resid <one> <two> <three> ...)"
#       Example: "protein and alpha and (resid 10 to 14)", "nucleic and name C1' and not (resid 11 15)"
set residue_selection 	"nucleic and name C1'"

# ======================= OUTPUT ===========================
set out_file_name 		"rmsf.dat";		# output file name
set out_delimiter 		" \t ";			# output delimiter 

set comment_token 		"#";	# Token used for Comments
set comment_header		1;		# Whether to comment out the columns header
# -----------------------------------------------------

puts "=============  RMSF  ================="

## Molecule ID (defaults to top, can be int)
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

### BASE SELECTION (for each frame) --------------
# Examples: "protein", "nucleic", "lipid", "all not water"

#puts "---------------- BASE SELECTION ---------------"
#puts " Used to select the part of system conating residues"
#puts " Examples: \"protein, \"nucleic\", \"lipid\", \"all not water\""
#puts ""
#puts -nonewline " -> set BASE SELECTION for each frame: "; flush stdout; 
#set base_selection [gets stdin]
#puts ""

#puts "--------------- RESIDUE SELECTION --------------"
#puts " Used to get RESIDUE_ID's of residues within BASE_SELECTION for which RMSF is required"
#puts " Examples: "
# FOR PROTEIN: "protein and alpha" -> get residue_id of each alpha C-atom in the protein (since there is only 1 alpha C per residue)
# FOR NUCLEIC Acids: "nucleic and name C1'" get residue_id of each C1'-atom of nucleotide sugar (since there is only 1 C1' per nucleotide)
#puts "   1. for Protein all residues: \"protein and alpha\" or \"protein and name CA\"" 
#puts "   2. for Nucleic all residues: \"nucleic and name C1'\""
#puts "   3. for RESIDUE RANGE, append \"and (resid <start> to <end>)\" or \"and (resid <one> <two>)\""
#puts "      Example: \"protein and alpha and (resid 10 to 14)\", \"nucleic and name C1'\" and not (resid 11 15)"
#puts ""
#puts -nonewline " -> set RESIDUE SELECTION: "; flush stdout; 
#set residue_selection [gets stdin]
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

puts "----------------- Frame Indices ------------------"
puts "INFO: Total Frames: $nf";
puts "Frame index must be in range \[0, ${last_index}\], or -ve for back indices"
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

log2file "${comment_token}================ RMSF =================="
puts "--------------------------"
log "${comment_token}LOG: INPUT Structure File: \"${psf_file}\" | Frame File(s): \[${frame_files}\]"
log "${comment_token}LOG: BASE SELECTION: \"${base_selection}\""
log "${comment_token}LOG: RESIDUE SELECTION: \"${residue_selection}\""
log "${comment_token}--------------------------"
log "${comment_token}LOG: Total Frames: ${nf}";
log "${comment_token}LOG: REFERENCE Frame Index: ${ref_frame_index}";
log "${comment_token}LOG: START Frame Index: ${start_frame_index} | END Frame Index: ${end_frame_index} | Frames for RMSF: [expr ${end_frame_index} - ${start_frame_index} + 1]"

set res_ids [[atomselect $mol_id $residue_selection] get resid];
set res_count [llength $res_ids];
log "${comment_token}LOG: Selected RESIDUE COUNT: $res_count";
puts "${comment_token}LOG: RESIDUE ID's: \[${res_ids}\]";
log "${comment_token}----------------------------"

# Reference Frame
set ref_frame [atomselect $mol_id $base_selection frame $ref_frame_index]

# Current frame (no frame yet)
set cur_frame [atomselect $mol_id $base_selection]

# make a selection with all atoms
set all_atoms [atomselect $mol_id all]

# Initialize rmsf array
foreach r $res_ids {
	set rmsf_arr($r) 0
}

# loop over all frames in the trajectory
for { set i $start_frame_index } { $i <= $end_frame_index } { incr i } {
	if {[expr $i % 10000] == 0} {
		puts "INFO: processing Frames $i-[expr min($i + 10000, $end_frame_index)]";
	}

	# get the current frame
	$cur_frame frame $i
	$all_atoms frame $i
	
	# compute the transformation
	set trans_mat [measure fit $cur_frame $ref_frame]
	# do the alignment
	$all_atoms move $trans_mat
	
	# compute the contribution to RMSF
	# loop through all residues
	foreach r $res_ids {
		set ref [atomselect $mol_id "${base_selection} and resid $r and noh" frame $ref_frame_index]
		set cur [atomselect $mol_id "${base_selection} and resid $r and noh" frame $i]
		set rmsf_arr($r) [expr $rmsf_arr($r) + [measure rmsd $cur $ref]]
		$cur delete;	# gc
		$ref delete;	# gc
	}
}

# Header for output file
if {$comment_header == 0} {
	set out_header "Frame${out_delimiter}RMSF";
} else {
	set out_header "${comment_token}Frame${out_delimiter}RMSF";
}

puts "";
log $out_header;

set average_rmsf 0
foreach r $res_ids {
    set rmsf_arr($r) [ expr $rmsf_arr($r) / $nf ]
    
    # print the RMSF and save in out_file
    log "$r${out_delimiter}$rmsf_arr($r)"
    
    # Computing average RMSF
    set average_rmsf [expr $average_rmsf + $rmsf_arr($r)]
    
    # Set the "user" field (of atoms of the residue) to RMSF value 
    set res_atoms [atomselect $mol_id "${base_selection} and resid $r"] 
    $res_atoms set user $rmsf_arr($r)
    $res_atoms delete;	# gc
}

set average_rmsf [ expr $average_rmsf / $res_count ]
puts ""
log "${comment_token}----------------------------"
log "${comment_token}INFO: Average RMSF per residue : ${average_rmsf}"

close $out_file

puts "=================  FINISHED  ===================="
puts "LOG: Output File: ${out_file_name}, delimiter: '$out_delimiter', comment token: '${comment_token}"
puts "================================================="

exit
