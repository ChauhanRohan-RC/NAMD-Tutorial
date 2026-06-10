#!/usr/bin/env -S vmd -dispdev text -e

##############################################################################################
### Script to extract PDB from binary Frame files (.coor, .dcd)             			       
###----------------------------
## USE CASES
# 1. BINARY FRAME (.coor) => PDB
# 2. Extract single FRAME from .dcd using FRAME_INDEX
# 3. SMD ATOM 1 -> 2 LINK VECTOR and SMD-PCF PDB file where transform_smd.tcl will not work (For ex. Periodic systems)
##############################################################################################

## USAGE ------------------------
# 1. Copy script to working dir
# 2. INPUT: Set input structure (.psf) and frame (.pdb, .coor) file
# 4. run with "./extract_pdb.tcl"
# 5. OUTPUT: generates .pdb file "<out_pdb_file>.pdb"


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


# INPUT ------------------
set psf_file			"../../common/dna_wb.psf";		# TODO: input strcuture file (.psf)
set frame_files		{ "dna_wb_eq.restart.coor" };		# TODO: input frame files (.pdb, .coor, .dcd)

# OUTPUT ------------------
set out_file_prefix 		"dna_wb_eq";				# TODO: [OPTIONAL] prefix for output files, "" for NONE

### BETA and OCCUPANCY -------
set set_beta		off;		# TODO: [on/off] whether to set beta values of all atoms
set beta_value	0;		# beta value to set (if $set_beta is ON)

set set_occupancy		off;		# TODO: [on/off] whether to set occupancy values of all atoms
set occupancy_value		0;		# occupancy value to set (if $set_occupancy is ON)

### LINK VECTOR b/w 2 atoms --------------
# -> Outputs Link Vector ATOM 1 -> ATOM 2
# -> Usually FIXED -> SMD atom
set calc_link_vector 	on;		# TODO: [on/off] Calculate Link vector between below atoms: Usually FIXED and SMD atoms
set selection_atom1		"nucleic and resid 1 and name C5'";		# selection for atom 1 (FIXED Atom)
set selection_atom2		"nucleic and resid 14 and name C3'";		# selection for atom 2 (SMD Atom)

set create_smd_pcf_file		on;		# only works when "calc_link_vector" = on
set smd_pcf_force			0.24;	# Force to write in smd-pcf file (in kcal/mol/A = pN / 69.5)




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


#
#
# =========================  MAIN  ==================================
# loading INPUT .psf and frame (.pdb, .coor) files
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

puts ""
puts ""
log "# ===============  EXTRACT PDB  =================="
log "# -> INPUT Structure File: \"${psf_file}\" | Frame File(s): \[${frame_files}\]"
log "# -> TOTAL Frames: $nf"
log "# -> SELECTED Frame Index: $frame_index"
log "# ----------------------------------------------------"
log "# -> INPUT Set BETA: $set_beta | BETA Value: $beta_value"
log "# -> INPUT Set OCCUPANCY: $set_occupancy | OCCUPANCY Value: $occupancy_value"
log "# ----------------------------------------------------"
puts ""

# Selecting Atoms
set alla [atomselect $mol_id "all" frame $frame_index];		# All AToms

# Setting BETA and OCCUPANCY --------------
puts ""
if {[string trim $set_beta] eq "on"} {
	puts "# -> LOG: SETTING BETA = $beta_value"
	$alla set beta $beta_value
}

if {[string trim $set_occupancy] eq "on"} {
	puts "# -> LOG: SETTING OCCCUPANCY = $occupancy_value"
	$alla set occupancy $occupancy_value
}
puts ""
# --------------------------------------

# Write output PDB files
set has_out_pdb 0
set has_out_smd_pcf_pdb 0

if { [info exists out_pdb_file] && [string trim $out_pdb_file] ne "" } {
	$alla writepdb $out_pdb_file
	
	set has_out_pdb 1
}


# Calculate Link Vector b/w 2 Atoms
if { [string trim $calc_link_vector] eq "on" && [string trim $selection_atom1] ne "" && [string trim $selection_atom2] ne "" } {
	set a1 [ atomselect $mol_id $selection_atom1 frame $frame_index ]
	set a2 [ atomselect $mol_id $selection_atom2 frame $frame_index ]
	
	set a1pos [lindex [$a1 get {x y z}] 0]
	set a2pos [lindex [$a2 get {x y z}] 0]
	
	set vec [vecsub $a2pos $a1pos];		# ATOM 1 -> ATOM 2 Link Vector
	set vec_unit [vecnorm $vec];			# ATOM 1 -> ATOM 2 Link UNIT Vector
	
	set dsq [vecdot $vec $vec];			# Distance Square (self dot)
	set dist [expr { sqrt($dsq) }];		# Distance in Å
	
	# SMD-PCF file
	if {[string trim $create_smd_pcf_file] eq "on" && [info exists out_smd_pcf_file] && [string trim $out_smd_pcf_file] ne "" && [info exists smd_pcf_force] } {
		# Normalize
		$alla set beta 0;
		$alla set occupancy 0;
		
		# FIXED ATOM
		$a1 set beta 1;
		
		# SMD ATOM
		$a2 set occupancy $smd_pcf_force;	# Force Magnitude
		$a2 set x [lindex $vec_unit 0];		# Force Direction
		$a2 set y [lindex $vec_unit 1];
		$a2 set z [lindex $vec_unit 2];
		
		# Write SMD-PCF file
		$alla writepdb $out_smd_pcf_file
		
		set has_out_smd_pcf_pdb	1
	}
	
	log ""
	log "### LINK VECTOR ------------------"
	log "# ATOM 1 (\"$selection_atom1\") => [$a1 get resname][$a1 get resid]:[$a1 get name]"
	log "# ATOM 2 (\"$selection_atom2\") => [$a2 get resname][$a2 get resid]:[$a2 get name]"
	log "# => ATOM 1 -> 2 DISTANCE : $dist Å"
	log "# => ATOM 1 -> 2 LINK UNIT VECTOR : ( $vec_unit )"
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

