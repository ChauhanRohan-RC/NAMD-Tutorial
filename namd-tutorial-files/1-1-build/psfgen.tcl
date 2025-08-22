############################################################################################
## Script to generate .psf (protein structure file) from a raw .pdb and topology file     ##
############################################################################################
#
# Outputs 2 files
#	1. new .pdb file: contains missing atoms (like H's) with guessed coordinates
#	2. .psf structure file
#
# USAGE
# 1. Copy this script in your working dir
# 2. set input and output config
# 2. In terminal, type "vmd -dispdev text -e psfgen.tcl"

## ------------- INPUT Config ---------------
set in_pdb "ubqp.pdb";		# TODO: input .pdb file
set topology_file "top_all27_prot_lipid.inp";	# TODO: TOPOLOGY file to use

## ------------- OUTPUT Config --------------
set out_name "ubq";			# TODO: Output .pdb and .psf file name
set out_segment_name "U";	# TODO: Segment name to use for ALL atoms of input .pdb file

## Main script
package require psfgen
topology $topology_file;

pdbalias residue HIS HSE;	# Alias for Histidine residue
pdbalias atom ILE CD1 CD;	# Alias for delta C atom of Isoleucine (ILE)

# creating a new segment in pdb format: segment <name> { <format> <in_file> } 
segment $out_segment_name {pdb $in_pdb};	
coordpdb $in_pdb $out_segment_name;		# loading .pdb coordinates into the segment

# Guess the coordinates of missing atoms (like H's that are usually not present in .pdb files)
# The guessed coords will be improved upon during energy-minimization phase
guesscoord;

# Writing output files
writepdb "$out_name.pdb";		# also contains missing atoms (like H's) with guessed coords
writepsf "$out_name.psf";
