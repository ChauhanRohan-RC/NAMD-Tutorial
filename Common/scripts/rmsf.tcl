## Script to calculate RMSF of each residue of a molecule over each timestep (TS) in VMD
## NOTE: the function also sets the "user" field of each atom to the RMSF value for its residue (Like a TAG)

### USAGE-------------------
# 1. Copy this script in working dir 
# 2. Load the .psf and .dcd in VMD
# 3. Go to Extensions > Tk Console
# 4. enter command: "source rmsf.tcl"
# 5. enter MOL_ID, BASE_SELECTION and RESIDUE_SELECTION
# 6. enter FRAME indices
# 7. Generates output file "rmsd.dat"

# ======================= CONFIG ===========================
set out_file_name "rmsf.dat";
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

log "#=============  RMSF  ================="

## Molecule ID (defaults to top, can be int)
puts -nonewline " -> set Molecule ID (default: top): "; flush stdout; 
set mol_id [gets stdin]
if { [string trim $mol_id] eq ""} { set mol_id top }
puts ""

### BASE SELECTION (for each frame) --------------
# Examples: "protein", "nucleic", "lipid", "all not water"

puts "---------------- BASE SELECTION ---------------"
puts " Used to select the part of system conating residues"
puts " Examples: \"protein, \"nucleic\", \"lipid\", \"all not water\""
puts ""
puts -nonewline " -> set BASE SELECTION for each frame: "; flush stdout; 
set base_selection [gets stdin]
puts ""

puts "--------------- RESIDUE SELECTION --------------"
puts " Used to get RESIDUE_ID's of residues within BASE_SELECTION for which RMSF is required"
puts " Examples: "
# FOR PROTEIN: "protein and alpha" -> get residue_id of each alpha C-atom in the protein (since there is only 1 alpha C per residue)
# FOR NUCLEIC Acids: "nucleic and name C1'" get residue_id of each C1'-atom of nucleotide sugar (since there is only 1 C1' per nucleotide)
puts "   1. for Protein all residues: \"protein and alpha\" or \"protein and name CA\"" 
puts "   2. for Nucleic all residues: \"nucleic and name C1'\""
puts "   3. for RESIDUE RANGE, append \"and (resid <start> to <end>)\" or \"and (resid <one> <two>)\""
puts "      Example: \"protein and alpha and (resid 10 to 14)\", \"nucleic and name C1'\" and not (resid 11 15)"
puts ""
puts -nonewline " -> set RESIDUE SELECTION: "; flush stdout; 
set residue_selection [gets stdin]
puts ""

### FRAMES INPUT ---------------------------
# Number of frames
set nf [molinfo $mol_id get numframes]
set last_index [expr $nf - 1]

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
puts "----------------------------"
log "#LOG: Mol ID: $mol_id";
log "#LOG: BASE SELECTION: '${base_selection}'"
log "#LOG: RESIDUE SELECTION: '${residue_selection}'"
log2file "#----------------------------"
log "#LOG: Frame Count: ${nf}";
log "#LOG: REFERENCE Frame Index: ${ref_frame_index}";
log "#LOG: START Frame Index: ${start_frame_index} | END Frame Index: ${end_frame_index} | Frames for RMSF: [expr ${end_frame_index} - ${start_frame_index} + 1]"

set res_ids [[atomselect $mol_id $residue_selection] get resid];
set res_count [llength $res_ids];
log "#LOG: Calculating RMSF for $res_count RESIDUES";
puts "LOG: RESIDUE ID's: $res_ids";
log "#----------------------------"

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
	if {[expr $i % 1000] == 0} {
		puts "INFO: processing Frames $i-[expr min($i + 1000, $end_frame_index)]";
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


# output header
set out_header "#RES_ID${out_delimiter}RMSF";
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
log "#----------------------------"
log "#INFO: Average RMSF per residue : $average_rmsf"

close $out_file
puts "";
puts "---------------------------------"
puts "FINISHED. Output File: ${out_file_name}, delimiter: '$out_delimiter'"
