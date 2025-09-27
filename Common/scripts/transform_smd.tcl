#!/usr/bin/env -S vmd -dispdev text -e

##############################################################################################
## Script to Transform a frame (from .pdb, .coor, .dcd) for SMD simulations
## ------------------------------------------------------------------------------------------
## NOTE: DO NOT USE in Periodic systems because cell dimensions depend on atomic coordinates					
## -> brings FIXED ATOM to origin															
## -> brings SMD ATOM to X-axis																
##############################################################################################

## USAGE ------------------------
# 1. Copy script to working dir
# 2. INPUT: Set input structure (.psf) and frame (.pdb, .coor, .dcd) file(s)
# 3. INPUT: Set selections for FIXED and SMD atom
# 4. run with "./transform_smd.tcl"
# 5. OUTPUT: generates .pdb file "<out_file_prefix>.pdb"
# 6. OUTPUT: if "create_smd_pcf_file = on", then generates "<out_file_prefix>.smd-pcf.pdb"

# ======================= INPUT ===========================
set psf_file			"../../common/dna.psf";				# TODO: input strcuture file (.psf)
set frame_files		{ "dna_gbis_eq.restart.coor" };		# TODO: input frame files (.pdb, .coor, .dcd)

set selection_atom_fixed		"nucleic and resid 1 and name C5'";	# TODO: selection for FIXED ATOM -> to be brought to origin
set selection_atom_smd		"nucleic and resid 14 and name C3'";		# TODO: selection for SMD ATOM -> to be brought to X-axis

# OUTPUT ------------------
set out_file_prefix 		"dna_gbis_eq";				# TODO: [OPTIONAL] prefix for output files, "" for NONE

### SMD File (OPTIONAL) --------------
# -> Calculates the Link Vector ATOM 1 (FIXED) -> ATOM 2 (SMD)
# -> Create a smd-pcf file with the link vector (x,y,z) and force (Occupancy) set on SMD atom 

set create_smd_pcf_file		on;		# only works when "calc_link_vector" = on
set smd_pcf_force			0.24;	# Force to write in smd-pcf file (in kcal/mol/A = pN / 69.5)

### BETA and OCCUPANCY -------
set set_beta		on;		# TODO: [on/off] whether to set beta values of all atoms
set beta_value	0;		# beta value to set (if $set_beta is ON)

set set_occupancy		on;		# TODO: [on/off] whether to set occupancy values of all atoms
set occupancy_value		0;		# occupancy value to set (if $set_occupancy is ON)





# OUTPUT FILES --------------------------
if { [info exists out_file_prefix] && [string trim $out_file_prefix] ne "" } {
	set out_pdb_file			$out_file_prefix.pdb;				
	set out_meta_file		$out_file_prefix.pdb.comments; # Store meta information
	set out_smd_pcf_file		$out_file_prefix.smd-pcf.pdb;
} else {
	set out_pdb_file			"";			
	set out_meta_file		"";
	set out_smd_pcf_file		"";
}

# =========================  MAIN  ==================================
# loading INPUT .psf and frame (.pdb, .coor, .dcd) files
set mol_id [mol new $psf_file];	

foreach frame_file $frame_files {
	puts "---------------------------------------------"
	puts " -> Loading Frames from file : ${frame_file}"
	puts "---------------------------------------------"
	mol addfile $frame_file waitfor all molid $mol_id;
}


### FRAMES INPUT ---------------------------
# Number of frames
set nf [molinfo $mol_id get numframes];
set last_index [expr $nf - 1];

if {$nf == 0} {
	puts "-------------------------------------------"
	puts "INFO: There are no frames in input file(s): \[${frame_files}\]"
	puts "LOG: QUITING (NO FRAMES)"
	puts "-------------------------------------------"
	exit
}

if {$nf == 1} {
	set frame_index 	0;
} else {
	puts "----------------- Frame Index ------------------";
	puts "INFO: Total Frames: $nf";
	puts "Frame index must be in range \[0, ${last_index}\], or -ve for back indices";
	puts ""
	puts -nonewline " -> Frame Index (default: 0): "; flush stdout; set frame_index [gets stdin]
	puts "------------------------------------------------"
	
	# NORMALIZING FRAME_INDEX (default 0)
	if { [string trim $frame_index] eq ""} {
		set frame_index 0
	} elseif {$frame_index < 0} {
		set frame_index [expr $last_index + (($frame_index + 1) % -$nf)]
	} elseif { $frame_index >= $nf } {
		set frame_index [expr $frame_index % $nf]
	}
} 


# META FILE -------
if { [info exists out_meta_file] && [string trim $out_meta_file] ne "" } {
	set meta_file [open $out_meta_file w];
	set has_meta 1
} else {
	set meta_file ""
	set has_meta 0
}

proc log2meta { msg } { 
	global has_meta;
	global meta_file;
	
	if {$has_meta == 1} {
		puts $meta_file $msg;	# to meta file
	}
}

proc log { msg } { 
	global has_meta;
	global meta_file;
	
	puts $msg; flush stdout;		# to stdout
	log2meta $msg;	# to meta file
}


# Selecting Atoms
set alla [atomselect $mol_id "all" frame $frame_index];		# All AToms
set fa [atomselect $mol_id $selection_atom_fixed frame $frame_index];		# Fixed-Atom
set sa [atomselect $mol_id $selection_atom_smd frame $frame_index];		# SMD-Atom

set fa_str "[$fa get resname][$fa get resid]:[$fa get name]"
set sa_str "[$sa get resname][$sa get resid]:[$sa get name]"

# LOGGING
if {1} {
	set fapos_pre [lindex [$fa get {x y z}] 0];
	set sapos_pre [lindex [$sa get {x y z}] 0];
	
	set vec_pre [vecsub $sapos_pre $fapos_pre];		# ATOM FIX -> ATOM SMD Link Vector
	set vec_unit_pre [vecnorm $vec_pre];			# ATOM FIX -> ATOM SMD Link UNIT Vector
	
	set dsq_pre [vecdot $vec_pre $vec_pre];		# Distance Square (self dot)
	set dist_pre [expr { sqrt($dsq_pre) }];		# Distance in Å

	puts ""
	puts ""
	log "# ===============   Transform for SMD  =================="
	log "# INPUT Structure File: \"${psf_file}\" | Frame File(s): \[${frame_files}\]"
	log "# INPUT FIXED ATOM selection: \"${selection_atom_fixed}\""
	log "# INPUT SMD ATOM selection: \"${selection_atom_smd}\""
	log "# ----------------------------------------------------"
	log "# INPUT Set BETA: $set_beta | BETA Value: $beta_value"
	log "# INPUT Set OCCUPANCY: $set_occupancy | OCCUPANCY Value: $occupancy_value"
	log "# ----------------------------------------------------"
	log "# -> TOTAL Frames: $nf"
	log "# -> SELECTED Frame Index: $frame_index"
	log "-------------------------"
	log "# -> FIXED Atom : ${fa_str} AT ( $fapos_pre )    (to be brought to origin)"
	log "# -> SMD Atom   : ${sa_str} AT ( $sapos_pre )    (to be brought to X-axis)"
	log ""
	log "### LINK VECTOR (Pre-Transform) ------------------"
	log "# => ATOM FIXED -> SMD  DISTANCE : $dist_pre Å"
	log "# => ATOM FIXED -> SMD  LINK UNIT VECTOR : ( $vec_unit_pre )"
	log ""

} 

# Bringing FIXED Atom to Origin 
set fapos [lindex [$fa get {x y z}] 0];		# FIXED atom position
$alla moveby [vecscale -1 $fapos];			# Bring FIXED atom to Origin

# Bringing SMD Atom to X-axis 
set sapos [lindex [$sa get {x y z}] 0];		# SMD atom position (after FIXED atom is brought to origin)
set trans_mat [transvecin $sapos];		# Transform Matrix to bring SMD atom to X-axis
$alla move $trans_mat;		# Apply transform
log "##########################"
log "# => TRANSFORM APPLIED !!!"
log "##########################"

#########################
# Write output file(s)
##########################

# Setting BETA and OCCUPANCY --------------
puts ""
if {[string trim $set_beta] eq "on"} {
	puts "# -> LOG: SETTING BETA = $beta_value"
	$alla set beta $beta_value;
}

if {[string trim $set_occupancy] eq "on"} {
	puts "# -> LOG: SETTING OCCCUPANCY = $occupancy_value"
	$alla set occupancy $occupancy_value;
}
puts ""
# --------------------------------------

# Write output PDB files
set has_out_pdb 0;
set has_out_smd_pcf_pdb 0;

if { [info exists out_pdb_file] && [string trim $out_pdb_file] ne "" } {
	$alla writepdb $out_pdb_file;
	
	set has_out_pdb 1;
}


# Calculate Link Vector b/w FIXED -> SMD ATOM and Create SMD-PCF file
if { 1 } {
	set fapos [lindex [$fa get {x y z}] 0];	# ATOM 1: FIXED
	set sapos [lindex [$sa get {x y z}] 0];	# ATOM 2: SMD
	
	set vec [vecsub $sapos $fapos];		# ATOM 1 -> ATOM 2 Link Vector
	set vec_unit [vecnorm $vec];			# ATOM 1 -> ATOM 2 Link UNIT Vector
	
	set dsq [vecdot $vec $vec];			# Distance Square (self dot)
	set dist [expr { sqrt($dsq) }];		# Distance in Å
	
	# SMD-PCF file
	if {[string trim $create_smd_pcf_file] eq "on" && [info exists out_smd_pcf_file] && [string trim $out_smd_pcf_file] ne "" && [info exists smd_pcf_force] } {
		# Normalize
		$alla set beta 0;
		$alla set occupancy 0;
		
		# FIXED ATOM
		$fa set beta 1;
		
		# SMD ATOM
		$sa set occupancy $smd_pcf_force;	# Force Magnitude
		$sa set x [lindex $vec_unit 0];		# Force Direction
		$sa set y [lindex $vec_unit 1];
		$sa set z [lindex $vec_unit 2];
		
		# Write SMD-PCF file
		$alla writepdb $out_smd_pcf_file;
		
		set has_out_smd_pcf_pdb	1;
	}
	
	log ""
	log "### LINK VECTOR (Post Transform) ------------------"
	log "# => ATOM FIXED -> SMD  DISTANCE : $dist Å"
	log "# => ATOM FIXED -> SMD  LINK UNIT VECTOR : ( $vec_unit )"
	puts "---------------------------------"
	puts ""
}


# LOGS --------------------
log ""
log "# -----------------------------------------"
if {$has_out_pdb == 1} {
	log "# -> OUTPUT PDB: \"${out_pdb_file}\""
}

if {$has_out_smd_pcf_pdb == 1} {
	log "# -> OUTPUT SMD-PCF PDB: \"${out_smd_pcf_file}\" | Force: $smd_pcf_force kcal/mol/Å"
}

log "# -----------------------------------------"

if { $has_meta == 1 } {
	puts "# -> OUTPUT META File: \"${out_meta_file}\""
	puts "# -----------------------------------------"
	puts ""
}

exit

