proc calc_ss_percentage {sel outfilename} {
    # Open the output file for writing
    set outfile [open $outfilename w]
    
    # Write the header to the file
    puts $outfile "Frame\tHelix(%)\tSheet(%)\tTurn(%)\tCoil(%)\tOther(%)"

    # Get the total number of frames in the top molecule
    set num_frames [molinfo top get numframes]

    # Select only the alpha-carbons (CA) of the protein. 
    # This ensures we only count 1 atom per residue.
    # set sel [atomselect top "protein and name CA"]
    set total_res [$sel num]

    # Check if protein exists
    if {$total_res == 0} {
        puts "Error: No protein CA atoms found."
        close $outfile
        return
    }

    puts "Starting calculation for $num_frames frames over $total_res residues..."

    # Loop through each frame
    for {set f 0} {$f < $num_frames} {incr f} {
        
        # Go to the current frame and update the display
        animate goto $f
        display update

        # Force VMD to recalculate secondary structure for this specific frame via STRIDE
        vmd_calculate_structure top

        # Update our atom selection for the current frame
        $sel frame $f
        $sel update

        # Get the secondary structure code for every residue
        set ss_list [$sel get structure]

        # Initialize counters
        set count_helix 0
        set count_sheet 0
        set count_turn  0
        set count_coil  0
        set count_other 0

        # Count the occurrences of each structure type
        # STRIDE Codes:
        # H (Alpha Helix), G (3_10 Helix), I (Pi Helix) -> Helix
        # E (Extended Sheet), B (Isolated Bridge) -> Sheet
        # T -> Turn
        # C -> Coil
        foreach ss $ss_list {
            if {$ss == "H" || $ss == "G" || $ss == "I"} {
                incr count_helix
            } elseif {$ss == "E" || $ss == "B"} {
                incr count_sheet
            } elseif {$ss == "T"} {
                incr count_turn
            } elseif {$ss == "C"} {
                incr count_coil
            } else {
                incr count_other
            }
        }

        # Calculate percentages (multiplied by 100.0 to enforce floating-point division)
        set p_helix [expr {($count_helix * 100.0) / $total_res}]
        set p_sheet [expr {($count_sheet * 100.0) / $total_res}]
        set p_turn  [expr {($count_turn * 100.0) / $total_res}]
        set p_coil  [expr {($count_coil * 100.0) / $total_res}]
        set p_other [expr {($count_other * 100.0) / $total_res}]

        # Format the numbers to 2 decimal places and write the row to the file
        set out_line [format "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f" $f $p_helix $p_sheet $p_turn $p_coil $p_other]
        puts $outfile $out_line
        
        # Print progress to the VMD console every 100 frames
        if {$f % 100 == 0} {
            puts "Processed frame $f..."
        }
    }

    # Clean up the atom selection and close the file
    $sel delete
    close $outfile
    puts "Done! Secondary structure data saved to: $outfilename"
}
