#!/usr/bin/env -S vmd -dispdev text -e

# ==========================================================
# Find residue sequence and unique residue names 
# output with pretty formatting
# ==========================================================

# -------------------------------
# INPUT
# -------------------------------

set psf_file 		"amyl_wb.psf";		# TODO
set frame_file		"amyl_wb.pdb";		# TODO

# TODO: base selection
set selection 		"protein";

set out_file 		"residues.txt";

# -------------------------------
# MAIN
# -------------------------------
set items_per_line 10;	# How many residues to print before wrapping to the next line
set col_width 6;		# Width of each column (6 ensures safe spacing)

set mol_id [mol new "$psf_file"];
mol addfile "$frame_file" first 0 last 0 molid $mol_id waitfor all;

# 1. Select the atoms you want to check (use "all" for the whole system)
set sel [atomselect $mol_id "${selection}"];

set res_indices [$sel get residue]
set res_names [$sel get resname]
set res_ids  [$sel get resid]

set seq_names {}
set seq_ids  {}
set current_res_idx -1

# Loop through all atoms and extract the name only when a new residue starts
foreach idx $res_indices name $res_names id $res_ids {
    if {$idx != $current_res_idx} {
        lappend seq_names $name
        lappend seq_ids $id
        set current_res_idx $idx
    }
}

set res_count [llength $seq_names];

# Sort the list and keep only the unique values
set unique_resnames [lsort -unique $seq_names];
set unique_res_count [llength $unique_resnames];

# -=-------------------------------
# OUTPUT 
# ---------------------------------
set out_fd [open "${out_file}" "w"];

proc log { msg } {
	global out_fd;
	
	puts $msg;
	puts $out_fd $msg;
}

# Print the result
log "\n\n"
log "INPUT PSF		 : $psf_file"
log "INPUT Frame     : $frame_file"
log "Selection       : $selection"
log "RESIDUE Count   : $res_count (total), $unique_res_count (unique)"
log "-----------------------------------------------------"
log "Unique Residues : $unique_resnames"
log "-----------------------------------------------------"

log "\nRESIDUE SEQUENCE\n"
# Loop through the lists in chunks
for {set i 0} {$i < $res_count} {incr i $items_per_line} {
    # Determine the end index of the current chunk
    set end [expr {$i + $items_per_line - 1}]
    if {$end >= $res_count} { 
        set end [expr {$res_count - 1}] 
    }
    
    # Grab the current slice of IDs and Names
    set chunk_ids   [lrange $seq_ids $i $end]
    set chunk_names [lrange $seq_names $i $end]
    
    set string_ids ""
    set string_names ""
    
    # Format each item to have exactly `col_width` spaces (left-aligned using %-Ns)
    foreach id $chunk_ids name $chunk_names {
        append string_ids   [format "%-${col_width}s" $id]
        append string_names [format "%-${col_width}s" $name]
    }
    
    # Print the stacked lines
    log $string_ids
    log $string_names
    log ""; # Add a blank line between chunks for readability
}
log "-----------------------------------------------------\n"



# Free up memory
$sel delete

close $out_fd;
exit;
