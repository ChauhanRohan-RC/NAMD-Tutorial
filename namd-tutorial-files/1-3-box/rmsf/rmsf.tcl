## Script to calculate RMSF of each residue of a molecule over each timestep (TS) in VMD
## NOTE: the function also sets the "user" field of each atom to the RMSF value for its residue (Like a TAG)

# 1. Load the .psf and .dcd in VMD
# 2. Go to Extensions > Tk Console
# 3. enter command: source rmsf.tcl

# 4. TO CALCULATE RMSF OF ALL RESIDUES: call function rmsf_all <mol> [out_file_name] [delimiter]
#--------------------- Example -----------------------------------------
#	source rmsf.tcl
#	rmsf_all 0 "rmsf_all.dat" "," (0: mol_id, optional file name and delimiter)
#	(generates output file "rmsf.dat" by default)
#-----------------------------------------------------------------------

# 5. TO CALCULATE RMSF OF SPECIFIC RESIDUES: call function rmsf <mol> <res_ids>
#--------------------- Example -----------------------------------------
#	source rmsf.tcl
#	set residue_ids {1 5 9}
#	rmsf top residue_ids rmsf_3.dat (optional file name and delimiter)
#	(generates output file "rmsf.dat" by default)
#-----------------------------------------------------------------------

# 5. TO CALCULATE RMSF OF SINGLE RESIDUE: call function rmsf_res <mol> <res_id>
#--------------------- Example -----------------------------------------
#	source rmsf.tcl
#	set residue_id 5
#	rmsf_res top residue_id
#-----------------------------------------------------------------------

proc all_res_ids {{mol top}} {
	# get residue_id of each alpha C-atom in the protein (since there is only 1 alpha C in each residue)
	set res_ids [[atomselect $mol "protein and alpha"] get resid];
	
	return $res_ids;
}

proc rmsf_all {{mol top} {out_file_name rmsf.dat} {delimiter " \t "}} {
	return [rmsf $mol [all_res_ids $mol] $out_file_name $delimiter]
}

###
# Calculates the RMSF of given residues over all frames
#
# @param mol: the molecule id, defaults to "top" if not given
# @param res_ids: array of residue ids
# @returns: array containing RMSF of the given residues over all simulation frames
###
proc rmsf {{mol top} res_ids {out_file_name rmsf.dat} {delimiter " \t "}} {
    #global out_file_name;
    
    # use frame 0 for the reference
    set ref_frame [atomselect $mol "protein" frame 0]
    # the current frame
    set cur_frame [atomselect $mol "protein"]
    # make a selection with all atoms
    set all_atoms [atomselect top all]
    
    # get the number of frames
    set nf [molinfo $mol get numframes]
    
    # Initialize rmsf array
    foreach r $res_ids {
		set rmsf_arr($r) 0
    }
    
    # loop over all frames in the trajectory
    for { set frame 0 } { $frame < $nf } { incr frame } {
    	puts "info: processing Frame $frame";
    
		# get the correct frame
		$cur_frame frame $frame
		$all_atoms frame $frame
		
		# compute the transformation
		set trans_mat [measure fit $cur_frame $ref_frame]
		# do the alignment
		$all_atoms move $trans_mat
		
		# compute the contribution to RMSF
		# loop through all residues
		foreach r $res_ids {
			set ref [atomselect $mol "protein and resid $r and noh" frame 0]
			set cur [atomselect $mol "protein and resid $r and noh" frame $frame]
			set rmsf_arr($r) [expr $rmsf_arr($r) + [measure rmsd $cur $ref]]
			$cur delete;	# gc
			$ref delete;	# gc
		}
    }
    
    
    # open output file for writing
    set out_file [ open $out_file_name w ]
    
    set average_rmsf 0
    
    puts "";
    puts "RES_ID${delimiter}RMSF";	# console header
	foreach r $res_ids {
	    set rmsf_arr($r) [ expr $rmsf_arr($r) / $nf ]
	    
	    # print the RMSF and save in out_file
	    set log "$r${delimiter}$rmsf_arr($r)"
	    puts "$log"
	    puts $out_file $log
	    
	    # Computing average RMSF
	    set average_rmsf [expr $average_rmsf + $rmsf_arr($r)]
	    
	    # Set the "user" field (of atoms of the residue) to RMSF value 
	    set res_atoms [atomselect $mol "resid $r"] 
        $res_atoms set user $rmsf_arr($r)
        $res_atoms delete;	# gc
	}
	
    set average_rmsf [ expr $average_rmsf / [llength $res_ids] ]
    puts "";
    puts "Average RMSF per residue : $average_rmsf"
    
    close $out_file
    puts "";
    puts "info: RMSF data saved to '$out_file_name' with delimiter '$delimiter'";
    
    #return [array get rmsf_arr]
}


###
# Calculates the RMSF of a given residue over all frames
# @param mol: the molecule id, defaults to "top" if not given
# @param res_id: residue id
# @returns: RMSF of a given residue over all simulation frames 
###
proc rmsf_res {{mol top} res_id} {
    
    # use frame 0 for the reference
    set ref [atomselect $mol "protein and resid $res_id and noh" frame 0]
    
    # the frame being compared
    set cur [atomselect $mol "protein and resid $res_id and noh"]

    #get the number of frames
    set nf [molinfo $mol get numframes]
    
    # RMSF
    set rmsf 0
    
    #loop over all frames in the trajectory
    for {set frame 1} {$frame < $nf} { incr frame } {
		# get the correct frame
		$cur frame $frame
		
		# compute the transformation and  do the alignment
		#$cur move [measure fit $cur $ref]
		
		# compute the RMSF contribution of this frame (RMSD of each atom within the residue)
		set rmsf [ expr $rmsf + [measure rmsd $cur $ref]]
    }
    
    set rmsf [ expr $rmsf / $nf ]
    puts "RMSF of residue $res_id : $rmsf"
}
