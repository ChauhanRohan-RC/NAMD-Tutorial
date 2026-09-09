#!/usr/bin/env -S vmd -dispdev text -e

# =======================================================================
# VMD TCL Script for Radial Distribution Function (RDF)
#-----------------------------------------------------------------
# => divides big dcd files in chunks and process sequentially
# NOTE: Much slower than AmberTools "cpptraj radial" command
# =======================================================================

# OUTPUT Columns
# 1. r        : radial distance from reference (Å) (bin center)
# 2. g_r      : radial distribution function RDF (normalized probability)
# 3. int_r    : running integral count (coordination number).
#               Avg number of selection_2 atoms within the sphere of radius r


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



#-----------------------------------
# INPUT Params (search TODO)
#-----------------------------------

# File Inputs
set psf_file		"../../common/amyl_wb.psf";		# TODO: input strcuture file (.psf)

# set dcd_files       { "run1.dcd" "run2.dcd" };    # TODO
set dcd_files  [find_files ".." "amyl_wb_eq" ".dcd" 3 3];  # <dir> <prefix> <suffix> [min_num] [max_num]

# Selections
# -----------------------------------
# 1. atom_selection1: reference subsystem (usually heavy atoms of protein/nucleic acid)
#    [protein heavy atoms] =>  "protein and noh"
#
# 2. atom_selection2: subsystem distributed around reference. Usually solvent (water Oxygen atoms)
#    [for water oxygen atom] => "water and name OH2"
#
set atom_selection1   "protein and noh";           # TODO:
set atom_selection2   "water and name OH2";        # TODO:

# RDF Parameters
set max_radius          10.0;        # Maximum distance for RDF calculation (Å)
set bin_size            0.1;         # Distance bin size (Å)
set use_pbc             1;           # Periodic boundary conditions. 1 = ON, 0 = OFF
set update_selections   0;           # update selections every frame [COSTLY]

# Performance & Memory Management
set chunk_size     1000;         # frames to process per chunk
set frame_step     1;          # skip frames during read


#------------------------
# OUTPUT
#------------------------
set out_file       "rdf.csv";
set out_delimiter  " ";
set comment_token  "#";



# =======================================================================
# 2. MAIN
# =======================================================================

puts "\n==================================================================="
puts "\[  INIT  \] Starting VMD Frame-Weighted Chunked RDF Calculator"
puts "\[ CONFIG \] Selection 1 : $atom_selection1"
puts "\[ CONFIG \] Selection 2 : $atom_selection2"
puts "\[ CONFIG \] Sel Update  : [expr {$update_selections ? "YES" : "NO"}] (every frame)"
puts "\[ CONFIG \] PBC Enabled : [expr {$use_pbc ? "YES" : "NO"}]"
puts "\[ CONFIG \] Max Radius  : $max_radius A | Bin Size: $bin_size A"
puts "\[ CONFIG \] Chunk Size  : $chunk_size | Frame Step: $frame_step"
puts "===================================================================\n"

set molid [mol new $psf_file waitfor all]

# Data arrays for frame-weighted accumulation
set num_bins          0;
set accumulated_gofr {};
set accumulated_int  {};
set r_values         {};
set total_chunks_processed 0;
set total_frames_processed 0;

set chunk_size [expr $chunk_size * $frame_step];

# Only for logging
set atom_count_sel1 0;
set atom_count_sel2 0;

# =======================================================================
# 3. CHUNKED PROCESSING LOOP
# =======================================================================

foreach dcd_file $dcd_files {
    puts "\[  READ  \] Opening trajectory: $dcd_file"
    set start_frame 0

    while {1} {
        set end_frame [expr $start_frame + $chunk_size - 1]

        # Load chunk into RAM, skipping frames dynamically to save I/O time
        animate read dcd $dcd_file beg $start_frame end $end_frame skip $frame_step waitfor all $molid

        # Determine exactly how many frames actually made it into RAM
        set num_loaded [expr [molinfo $molid get numframes] - 1]
        if {$num_loaded <= 0} { break }

        set actual_end [expr $start_frame + ($num_loaded * $frame_step) - 1]
        puts "\n\[ PROGRESS \] -> Processing raw frames $start_frame - $actual_end (Loaded $num_loaded frames)...\n"

        # Define selections for this chunk
        set sel1 [atomselect $molid $atom_selection1]
        set sel2 [atomselect $molid $atom_selection2]

        # Calculate RDF. We use step 1 here because skipping was already handled by 'animate read'
        set rdf_data [measure gofr $sel1 $sel2 delta $bin_size rmax $max_radius usepbc $use_pbc selupdate $update_selections step 1]

        # Initialization Logic
        if {$total_chunks_processed == 0} {
            set r_values [lindex $rdf_data 0];

            set num_bins [llength $r_values];
            for {set i 0} {$i < $num_bins} {incr i} {
                lappend accumulated_gofr 0.0;
                lappend accumulated_int  0.0;
            }

            set atom_count_sel1 [$sel1 num];
            set atom_count_sel2 [$sel2 num];

            puts "\n----------------------------------------------------"
            puts "Selection 1: \"$atom_selection1\"  ($atom_count_sel1 atoms)"
            puts "Selection 2: \"$atom_selection2\"  ($atom_count_sel2 atoms)"
            puts "----------------------------------------------------\n"
        }

        # Extract g(r) and the Integral of g(r)
        set chunk_gofr [lindex $rdf_data 1]
        set chunk_int  [lindex $rdf_data 2]

        # Frame-Weighted Accumulation (mathematically prevents skewed averages from short tail-chunks)
        set temp_g {}
        set temp_i {}
        foreach cur_g $accumulated_gofr new_g $chunk_gofr cur_i $accumulated_int new_i $chunk_int {
            lappend temp_g [expr {$cur_g + ($new_g * $num_loaded)}]
            lappend temp_i [expr {$cur_i + ($new_i * $num_loaded)}]
        }
        set accumulated_gofr $temp_g
        set accumulated_int  $temp_i

        # Prevent memory leaks
        $sel1 delete
        $sel2 delete
        animate delete beg 1 end -1 $molid

        incr total_chunks_processed
        incr total_frames_processed $num_loaded
        set start_frame [expr $start_frame + $chunk_size]
    }
}

# =======================================================================
# 4. DATA EXPORT & METADATA WRITING
# =======================================================================

if {$total_frames_processed > 0} {
    puts "\n\[  CALC  \] Averaging over $total_frames_processed analyzed frames..."
    set out [open $out_file w]

    # Write metadata block at the top
    puts $out "${comment_token} ==================================================================="
    puts $out "${comment_token} VMD Radial Distribution Function (RDF) & Integral Calculation"
    puts $out "${comment_token} ==================================================================="
    puts $out "${comment_token} Structure File : $psf_file"
    puts $out "${comment_token} Trajectories   : $dcd_files"
    puts $out "${comment_token} Selection 1    : \"$atom_selection1\"  ($atom_count_sel1 atoms, reference subsystem)"
    puts $out "${comment_token} Selection 2    : \"$atom_selection2\"  ($atom_count_sel2 atoms, distributed subsystem)"
    puts $out "${comment_token} Update Sel     : [expr {$update_selections ? "Enabled" : "Disabled"}] (every frame)"
    puts $out "${comment_token} ------------------------------------------------------------------"
    puts $out "${comment_token} Max Radius     : $max_radius Å"
    puts $out "${comment_token} Bin Size       : $bin_size Å"
    puts $out "${comment_token} Periodic Bound : [expr {$use_pbc ? "Enabled" : "Disabled"}]"
    puts $out "${comment_token} ------------------------------------------------------------------"
    puts $out "${comment_token} Frame Step     : $frame_step"
    puts $out "${comment_token} Total frames   : $total_frames_processed (processed)"
    puts $out "${comment_token} ------------------------------------------------------------------"
    puts $out "${comment_token} OUTPUT Columns"
    puts $out "${comment_token} => r        : radial distance from reference (Å) (bin center)"
    puts $out "${comment_token} => g_r      : radial distribution function RDF (normalized probability)"
    puts $out "${comment_token} => int_r    : running integral count (coordination number). Avg number of selection_2 atoms within the sphere of radius r"
    puts $out "${comment_token} ==================================================================="

    # Write CSV header
    puts $out "r${out_delimiter}g_r${out_delimiter}int_r"

    # Write data rows
    foreach r $r_values g_sum $accumulated_gofr i_sum $accumulated_int {
        # Divide the weighted sums by the total number of frames processed
        set g_avg [format "%.4f" [expr {$g_sum / double($total_frames_processed)}]]
        set i_avg [format "%.4f" [expr {$i_sum / double($total_frames_processed)}]]
        set r_fmt [format "%.4f" $r]

        puts $out "$r_fmt${out_delimiter}$g_avg${out_delimiter}$i_avg"
    }

    close $out
    puts "==================================================================="
    puts "\[ SUCCESS \] Calculation complete!"
    puts "\[ SUCCESS \] Result saved to: >>>  $out_file  <<<"
    puts "===================================================================\n"
} else {
    puts "\[ ERROR \] No frames were successfully processed. Check your DCD paths or frame step."
}

quit
