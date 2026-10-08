#!/usr/bin/env -S vmd -dispdev text -e

# =========================================================================================
# Finds the closest Water molecule to a given reference selection in a SINGLE FRAME
# =========================================================================================


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


# Read atomselection string from file (must contain only a single line)
proc read_selection_from_file { file_path } {
	set fd [open $file_path];
	set data [read $fd];
	close $fd;

	set trimmed [string trim $data];
	return $trimmed;
}


# -----------------------
# INPUT
# -----------------------
set psf_file		"../../common/amyl_wb.psf";					# TODO: input strcuture file (.psf)

# TODO: Single Frame File (.pdb, .coor) or DCD (only first frame is read)
set frame_file	[lindex [find_files "../" "amyl_wb_eq" ".dcd"] 0];  # <dir> <prefix> <suffix> [min_num] [max_num]

set reference_selection		[read_selection_from_file "selection_prot.txt"];		# TODO: find closest distance from this
set water_oxygens			"water and noh";									
set cutoff         			6.0;				# Search radius in Angstroms to optimize speed

# -----------------------
# OUTPUT
# -----------------------
set out_selection_file		"selection_water.txt"


# ======================================
# MAIN
# ======================================
# 0. Load molecule and frame
set mol_id [mol new "$psf_file" waitfor all];
mol addfile $frame_file first 0 last 0 waitfor all molid $mol_id;

# 1. Select the target protein atom
set target [atomselect $mol_id $reference_selection]
if {[$target num] != 1} {
    puts "Error: Protein selection must result in exactly 1 atom. Found [$target num] atoms."
    $target delete
    return
}
set target_coord [lindex [$target get {x y z}] 0]

# 2. Select nearby water oxygens to act as search candidates
set candidates [atomselect $mol_id "($water_oxygens) and within $cutoff of ($reference_selection)"]
if {[$candidates num] == 0} {
    puts "Error: No waters found within ${cutoff}Å."
    $target delete
    $candidates delete
    return
}

# 3. Iterate through candidates to find the minimum distance
set min_dist 99999.0
set closest_idx -1

foreach idx [$candidates get index] coord [$candidates get {x y z}] {
    set dist [veclength [vecsub $target_coord $coord]]
    if {$dist < $min_dist} {
        set min_dist $dist
        set closest_idx $idx
    }
}

set closest_atom [atomselect $mol_id "index $closest_idx"]
set c_resid   [lindex [$closest_atom get resid] 0]
set c_segname [lindex [$closest_atom get segname] 0]
set c_chain   [lindex [$closest_atom get chain] 0]

# 4. Build the robust selection string dynamically
set final_water_sel "water and resid $c_resid"

# Add segname if it exists in your topology (common in CHARMM/NAMD)
if {$c_segname != ""} {
    append final_water_sel " and segname $c_segname"
}
# Add chain if it exists in your topology (common in PDBs)
if {$c_chain != ""} {
    append final_water_sel " and chain $c_chain"
}

set min_dist_str [format "%.2f" $min_dist];

# Free up memory
$target delete
$candidates delete

# Output file
set fh [open "${out_selection_file}" w];
puts -nonewline $fh "${final_water_sel}";
close $fh

puts "---------------------------------------------------"
puts "Closest water distance : $min_dist_str Å"
puts "OUTPUT Selection       : \"$final_water_sel\""
puts "OUTPUT Selection File  : $out_selection_file"
puts "---------------------------------------------------"

exit;
