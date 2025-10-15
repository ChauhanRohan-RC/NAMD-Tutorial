#!/usr/bin/env -S vmd -dispdev text -e

############################################################################################
## Script to generate .psf (protein structure file) from a raw .pdb and topology file     ##
############################################################################################

# OUTPUT: 2 files
#	1. new .pdb file: contains missing atoms (like H's) with guessed coordinates
#	2. .psf structure file

# USAGE
# 1. Copy this script in your working dir
# 2. set input and output config
# 2. run with "./psfgen.tcl"


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


# ======================
# INPUT
# ======================
set input_pdb "ubqp.pdb";		# TODO: input .pdb file
set topology_files { "top_all27_prot_lipid.inp" };	# TODO: TOPOLOGY files to use

# ======================
# OUTPUT
# ======================
set out_name "ubq";			# TODO: Output .pdb and .psf file name
set out_segment_name "U";	# TODO: Segment name to use for ALL atoms of input .pdb file



# ======================
# MAIN
# ======================
package require psfgen;
foreach top_file $topology_files {
	topology $top_file;
}

# Aliases (Compatibility with AMBER and others)
pdbalias residue HID HSD;	# Alias for Histidine residue (Neutral, H on delta N)
pdbalias residue HIE HSE;	# Alias for Histidine residue (Neutral, H on epsilon N)
pdbalias residue HIP HSP;	# Alias for Histidine residue (Protonated, proton on epsilon N)
pdbalias atom ILE CD1 CD;	# Alias for delta C atom of Isoleucine (ILE)

# creating a new segment in pdb format: segment <name> { <format> <in_file> } 
segment $out_segment_name {pdb $input_pdb};	
coordpdb $input_pdb $out_segment_name;		# loading .pdb coordinates into the segment

# Guess the coordinates of missing atoms (like H's that are usually not present in .pdb files)
# The guessed coords will be improved upon during energy-minimization phase
guesscoord;

# Writing output files
writepdb "$out_name.pdb";		# also contains missing atoms (like H's) with guessed coords
writepsf "$out_name.psf";

exit;
