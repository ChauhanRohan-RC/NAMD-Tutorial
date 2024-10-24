##############################################################################################
### Script to Transform a frame file (.pdb, .coor) for SMD simulations		   				##
## -> brings FIXED ATOM to origin															##
## -> brings SMD ATOM to X-axis																##
##############################################################################################

## USAGE ------------------------
# 1. Copy script to working dir
# 2. INPUT: Set input structure (.psf) and frame (.pdb, .coor) file
# 3. INPUT: Set selections for FIXED and SMD atom
# 4. run with "vmd -dispdev text -e transform_smd.tcl"
# 5. OUTPUT: generates .pdb file "<out_pdb_file>.pdb"

# ======================= INPUT ===========================
set psf_file		"dna.psf";							# TODO: input strcuture file (.psf)
set in_frame_file	"dna_gbis_eq4.restart.coor";		# TODO: input frame file (.pdb, .coor)

set selection_atom_fixed	"resid 16 and name C3'";	# TODO: selection for FIXED ATOM -> to be brought to origin
set selection_atom_smd		"resid 1 and name C5'";		# TODO: selection for SMD ATOM -> to be brought to X-axis

# ======================= OUTPUT ===========================
set out_pdb_file 		"dna_gbis_eq4.pdb";				# TODO: output .pdb file

# =========================  MAIN  ==================================
# loading INPUT .psf and .dcd files
set mol_id [mol new $psf_file];	
mol addfile $in_frame_file waitfor all molid $mol_id;

# Selecting Atoms
set alla [atomselect $mol_id "all"];		# All AToms
set fa [atomselect $mol_id $selection_atom_fixed];		# Fixed-Atom
set sa [atomselect $mol_id $selection_atom_smd];		# SMD-Atom

set fa_str "[$fa get resname][$fa get resid]:[$fa get name]"
set sa_str "[$sa get resname][$sa get resid]:[$sa get name]"

puts "==============  Transform for SMD  ================="
puts "# INPUT frame file: \"${in_frame_file}\""
puts "# INPUT FIXED ATOM selection: \"${selection_atom_fixed}\""
puts "# INPUT SMD ATOM selection: \"${selection_atom_smd}\""
puts "-------------------------"
puts " -> FIXED Atom: ${fa_str}   (to be brought to origin)"
puts " -> SMD Atom: ${sa_str}     (to be brought to X-axis)"
puts "-------------------------"

# Bringing FIXED Atom to Origin 
set fapos [lindex [$fa get {x y z}] 0];		# FIXED atom position
$alla moveby [vecscale -1 $fapos];			# Bring FIXED atom to X-axis

# Bringing SMD Atom to X-axis 
set sapos [lindex [$sa get {x y z}] 0];		# SMD atom position (after FIXED atom is brought to origin)
set trans_mat [transvecin $sapos];		# Transform Matrix to bring SMD atom to X-axis
$alla move $trans_mat;		# Apply transform

# Write output file
$alla writepdb $out_pdb_file
flush stdout;

puts "=================  FINISHED  ===================="
puts "LOG: Output File: \"${out_pdb_file}\""
puts "================================================="

exit

