#!/usr/bin/env -S vmd -dispdev text -e

#################################################################################################
### Script to calculate Energies from simulation trajectories (DYNAMIC SELECTIONS)			   ##
# ---------------------------------------------------------------------------------------------##
# This script is best for DYNAMIC atom selections (which update every frame)
# For STATIC SELECTIONS, use namd_energy.static.sh (much faster)
#----------------------------------------------------------------------------------------------##
# Custom implementation of NAMD-Energy plugin of VMD to make it effecient and fast			   ##
# Created by RC, Feb 16 2026																   ##
#----------------------------------------------------------------------------------------------##
# -> Calculates Energies (Forces) on a subset of the system (due to itself/any other parts)    ##
# -> Calculate self and cross interaction energies of a subset of the system			       ##
# ---------------------------------------------------------------------------------------------##
# This script is highly efficient than NAMD Energy Plugin of VMD.
# It does not load dcd files into memory, and processes trajectory frame-by-frame.
#
## LIMITATIONS -----------
# 1. Only calculate the absolute force between two selections.
#	 CANNOT calculate the force projection on the vector connection COM of two selections, which
#	 is anyways not required in most cases
#################################################################################################


## USAGE --------------------------------------------------
# 0: First run normal simulation to obtain .dcd trajectories
# 1. Copy script to working dir
# 2. INPUT: Set input strcuture (.psf) and trajectories (.dcd)
# 3. INPUT: Set selection 1, Selection 2 (Optional), out_energies (Optional), out_file_prefix
#-----------------------NAMD Config file---------------------
# -> ALTERNATIVELY, SET ENVIRONMENT VARIABLES (used when variables are not set in script)
#--------------------------------------------
#	-> NAMD_ENERGY_SELECTION1	=	selection1
#	-> NAMD_ENERGY_SELECTION2	=	selection2 		  (optional)
#	-> NAMD_ENERGY_OUT_ENERGIES	= 	out_energies	  (optional)
#	-> NAMD_ENERGY_OUT_PREFIX	=	out_file_prefix
#	-> NAMD_ENERGY_TIMESTEP_END	=	timestep_end      (optional)
#	-> NAMD_ENERGY_LABEL	    =	label			  (optional)
#   -> NAMD_ENERGY_PROCESSES    =   namd_processes    (optional)
#--------------------------------------------
# 4. set other input and output params [search for TODO]
# 5. run with "./namd_energy.dynamic.tcl"
# 	OR
# 6. use namd_energy.sh launcher.
#	-> First, unset selection1, selection2, out_energies, out_file_prefix in this script
#	-> Set environment variables in namd_energy.sh
#		=> ./namd_energy.sh

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

# Total number of threads to run on
set namd_processes  8;
catch { set namd_processes "$::env(NAMD_ENERGY_PROCESSES)" }; # from env-variable


# NAMD Command
set namd_pes 		1;	# threads per NAMD process
set namd_cmd_list	[list "$::env(NAMD_MULTICORE)/namd3" "+p${namd_pes}"];	# more threads does not do anything
set debug			0;	# debug mode

# =============================
# INPUT
# =============================
# TODO: [OPTIONAL] Force-Field Parameter files
set param_files	{ "../../common/ff/par_all36m_prot.prm" "../../common/ff/toppar_water_ions.prot.str" };
set drude	off;		# Drude additive force field

# TODO: input strcuture file (.psf)
set psf_file		"../../common/amyld_wb.psf";

# TODO: LIST of frames (.dcd files) separated by space
# set dcd_files	{ "../amyl_wb_eq.dcd" };
set dcd_files	[find_files ".." "amyld_wb_eq" ".dcd"];  # <dir> <prefix> <suffix> [min_num] [max_num]

# Selection
set selection1	"";			# or set by  ENV VAR: NAMD_ENERGY_SELECTION1
set selection2	"";			# or set by  ENV VAR: NAMD_ENERGY_SELECTION2

## Output file names
set out_file_prefix 		"";		# or set by  ENV VAR: NAMD_ENERGY_OUT_PREFIX

## [OPTIONAL][ Label for this run
set label_		"";			# or ser by  ENV VAR: NAMD_ENERGY_LABEL

### Energies to calculate (as sequence of 4-letter codes)
# -------------------------------------------------------------------------------------
# OPTIONS: -bond -angl -dihe -impr -conf -vdw -elec -nonb -boun -misc -kine -pote -all
# -------------------------------------------------------------------------------------
# -> -boun (boundary), -kine (kinetic), -pote (potential)
# -> -conf (confomational) = bond + angle + dihedral + improper
# -> -nonb (non-nonded)    = elec + vdw
# -> -pote (potential)     = bond + angl + dihe + impr + vdw + elec + boun + misc
# -------------------------------------------------------------------------------------
# WITH SELECTION 2: ONLY [ -vdw -elec -nonb -boun -misc -pote -all ] ARE ALLOWED

set out_energies			"";		# or set by  ENV VAR: NAMD_ENERGY_OUT_ENERGIES
								# example: "-bond -elec -kine - pote"


## Params
set temperature			300;		# Temperature (in K) used for main simulation [-T]
set cutoff				12;		# Cutoff distance (in Å)		[-cutoff]
set switchdist			10;		# Switch distance (in Å) for non-bonded interactions [-switch]
								#   Normally < cutoff. -ve to turn off switching
set dielectric			1.0;		# Dielectric constant (> 1 will lessen the electrostatic forces)

## Periodic [OPTIONAL]
set initial_ext_sys		"../../wb_eq/amyld_wb_eq51.restart.xsc";		# [ONLY PERIODIC] Initial cell dimensions in .dcd file
set pme					on;		# [on/off] [ONLY PERIODIC] PME for long-range electrostatics. Only work when <initial_ext_sys> is specified
set pme_grid_spacing		1.0;		# (in Å) spacing b/w PME grid points on cell basis vectors, used for automatic PME grid sizes


## TODO: TIme Step parameters (ONLY USED FOR OUTPUT COLUMNS, DOES NOT AFFECT CALCULATION)
set timestep_first		0;			# First timestep.						[-ts]
set frame_freq			100;		# Timesteps b/w each frame = dcdfreq	[-stride]

# Timstep to stop the calculation at [INCLUSIVE, -1 for None]
set timestep_end 		-1;		# or set by  ENV VAR: NAMD_ENERGY_TIMESTEP_END

# ==================================
# OUTPUT Params
# ==================================

## Force Output	(ONLY APPLICABLE when SEL-2 is defined)
# -> Calculates force on SEL-1 due to SEL-2
set out_force			on;		# [on/off] show force b/w two selections. [-keepforce]
# set out_force_proj		off;		# [on/off] Only show projection of force on the vector connecting COM of two selections. [-projforce].  [output force will be signed (+ve: replusive, -ve attractive)]

set frame_skip			0;		# Calculate energy every <frame_skip> frame [-skip]

set out_delimiter 		" ";
set out_energy_format 	"%.4f";
set comment_token 		"#";			# For Comments. "" to disable comments




# =========================================================
# MAIN
# =========================================================

set time_start [clock seconds];

#------------------
# Init and basic checks
#------------------
if {[llength $dcd_files] == 0} {
	puts "\n==================================================="
    puts "-> ERROR: No dcd files specified !!!"
    puts "===================================================\n"
	exit;
}

# Output file prefix
if { [info exists out_file_prefix] == 0 || [string trim $out_file_prefix] eq "" } {

	# Trying to find env variable
	set failed [catch { set out_file_prefix "$::env(NAMD_ENERGY_OUT_PREFIX)" }];
	if { $failed != 0 } {
		# Environment variable not found
	}
	unset failed;


	if { [info exists out_file_prefix] == 0 || [string trim $out_file_prefix] eq "" } {
		puts "\n==================================================="
		puts "-> ERROR: Output file prefix must be defined !!!"
		puts "===================================================\n"
		exit;
    }
}

if { [info exists selection1] == 0 || [string trim "$selection1"] eq "" } {

	# Trying to find env variable
	set failed [catch { set selection1 "$::env(NAMD_ENERGY_SELECTION1)" }];
	if { $failed != 0 } {
		# Environment variable not found
	}
	unset failed;

	if { [info exists selection1] == 0 || [string trim "$selection1"] eq "" } {
		puts "\n==================================================="
		puts "-> ERROR: SELECTION 1 must be defined !!!"
		puts "===================================================\n"
		exit;
	}
}

# Selection 2-------------
set has_sel2 0;
if { [info exists selection2] == 0 || [string trim "$selection2"] eq "" } {

	# Trying to find env variable
	set failed [catch { set selection2 "$::env(NAMD_ENERGY_SELECTION2)" }];
	if { $failed != 0 } {
		# Environment variable not found
	}
	unset failed;
}

if { [info exists selection2] && [string trim "$selection2"] ne "" } {
	set has_sel2 1;

	if { [string trim "$selection1"] eq [string trim "$selection2"] } {
        puts "\n==================================================="
        puts "-> ERROR: SELECTION 1 and 2 MUST BE DIFFERENT (given \"$selection1\") !!!"
        puts "===================================================\n"
        exit;
	}
}

# Output energies parsing -----------------

# Try env variable
if { [info exists out_energies] == 0 || [string trim $out_energies] eq "" } {
    set failed [catch { set out_energies "$::env(NAMD_ENERGY_OUT_ENERGIES)" }];
	if { $failed != 0 } {
		# Environment variable not found
	}
	unset failed;
}

# If out_energies not specified, set to -all
if { [info exists out_energies] == 0 || [string trim $out_energies] eq "" || [string trim $out_energies] eq "-all"} {
    if { $has_sel2 == 1 } {
        set out_energies "-nonb, -boun, -misc, -pote";		# Only non-bonded allowed b/w two selections
    } else {
        set out_energies "-all";
    }
}

# Only Non bonding energies are allowed b/w 2 selections
if { $has_sel2 == 1 } {
	if { [regexp {bond|angl|dihe|impr|conf|kine} $out_energies] } {
		puts "\n=============================================================="
		puts "-> ERROR: Only Non-Bonded energies (-vdw, -elec, -nonb, -boun, -misc, -pote) can be calculate for interaction between two selections (given \"$out_energies\") !!!"
		puts "==============================================================\n"
		exit;
	}
}


# Perioidc BC
set has_ext_sys 0;
if { [info exists initial_ext_sys] && [string trim $initial_ext_sys] ne "" } {
	set has_ext_sys 1;
}

set has_pme 0;
if { [info exists pme] && [string trim $pme] eq "on" } {
	set has_pme 1;
}


if { [info exists frame_skip] && [string trim $frame_skip] ne "" && $frame_skip > 0 } {
	set has_skip 1;
} else {
	set has_skip 0;
	set frame_skip 0;
}

# Only when selection 2 is also defined
set keepforce 0;
set keepforce_str "OFF"
if { $has_sel2 == 1 && [info exists out_force] && [string trim $out_force] eq "on" } {
	set keepforce 1;
	set keepforce_str "ON";
}

set has_drude 0;
if { [info exists drude] && [string trim $drude] eq "on" } {
	set has_drude 1;
}

# set projforce 0;
# if { $keepforce == 1 && [info exists out_force_proj] && [string trim $out_force_proj] eq "on" } {
# 	set projforce 1;
# }


# Timestep End
if { [info exists timestep_end] == 0 || [string trim $timestep_end] eq "" || $timestep_end == -1 } {
    set failed [catch { set timestep_end "$::env(NAMD_ENERGY_TIMESTEP_END)" }];
	if { $failed != 0 } {
		# Environment variable not found
		set timestep_end -1;
	}

	unset failed;
}

## Label
if { [info exists label_] == 0 || [string trim $label_] eq "" } {
	# Trying to find env variable
	set failed [catch { set label_ "$::env(NAMD_ENERGY_LABEL)" }];
	if { $failed != 0 } {
		# Environment variable not found
	}
	unset failed;

	if { [info exists label_] == 0 || [string trim $label_] eq "" } {
		# Deafult Label
		if { $has_sel2 == 1 } {
			if { $keepforce == 1 } {
				set label_ "Interaction Energy and Forces (on Selection-1 due to Selection-2)";
			} else {
				set label_ "Interaction Energy b/w Selection-1 and Selection-2";
			}
		} else {
			set label_ "Self-Interaction Energy of Selection-1";
		}
	}
}




#-----------------------------------
# Atom Checks
#-----------------------------------
set mol_id [mol new $psf_file];
# Load the first frame as template
mol addfile [lindex $dcd_files 0] type dcd first 0 last 0 waitfor all molid $mol_id;

set all [atomselect $mol_id "all"];
set atom_count_total [$all num];
if {$atom_count_total == 0} {
	puts "\n==================================================="
	puts "-> ERROR: No Atoms found in DCD file: [lindex $dcd_files 0] !!!"
	puts "===================================================\n"
	exit;
}

set sel1 [atomselect $mol_id "$selection1"];
set atom_count_sel1 [$sel1 num];
if { $atom_count_sel1 == 0 } {
	puts "\n============================================================"
	puts "-> WARNING: No atoms found for selection 1 in first frame: \"$selection1\" !!!"
	puts "============================================================\n"
# 	exit;
}


if { $has_sel2 == 1 } {
	set sel2 [atomselect $mol_id "$selection2"];
	set atom_count_sel2 [$sel2 num];

	if { $atom_count_sel2 == 0 } {
		puts "\n============================================================"
		puts "-> WARNING: No atoms found for selection 2 in first frame: \"$selection2\" !!!"
		puts "============================================================\n"
# 		exit;
	}
} else {
	set sel2 "";
	set atom_count_sel2 0;
}



#---------------------------------------
# Worker Pool
#---------------------------------------
namespace eval workerpool {
    variable max_workers 3
    variable running 0
    variable next_id 0
    variable queue {}
    variable jobs
    variable idle_callback {}
    variable done 0
    variable use_shell 0

    # --------------------------
    # Init
    # --------------------------
    proc init {size {shell 0}} {
        variable max_workers
        variable use_shell
        variable running
        variable queue
        variable jobs
        variable done

        set max_workers $size
        set use_shell $shell
        set running 0
        set queue {}
        array unset jobs
        set done 0
    }

    proc on_idle {callback} {
        variable idle_callback
        set idle_callback $callback
    }

    proc capacity {} {
        variable running
        variable max_workers
        return [expr {$max_workers - $running}]
    }

    # --------------------------
    # Submit (WITH OPTIONAL TAG)
    # --------------------------
    proc submit {cmd callback {tag ""}} {
        variable next_id
        variable queue

        incr next_id
        set id $next_id

        lappend queue [list $id $cmd $callback $tag]
        run_next
        return $id
    }

    # --------------------------
    # Dispatcher
    # --------------------------
    proc run_next {} {
        variable running
        variable max_workers
        variable queue

        while {$running < $max_workers && [llength $queue] > 0} {
            set job [lindex $queue 0]
            set queue [lrange $queue 1 end]
            lassign $job id cmd callback tag
            start_job $id $cmd $callback $tag
        }

        check_done
    }

    # --------------------------
    # Start Process
    # --------------------------
    proc start_job {id cmd callback tag} {
        variable running
        variable jobs
        variable use_shell

        incr running

        if {$use_shell} {
            if {$::tcl_platform(platform) eq "windows"} {
                set pipe [open [list | cmd /c $cmd] r]
            } else {
                set pipe [open [list | /bin/sh -c $cmd] r]
            }
        } else {
            set pipe [open [linsert $cmd 0 |] r]
        }

        fconfigure $pipe -blocking 0 -buffering none

        set jobs($id,pipe) $pipe
        set jobs($id,cmd) $cmd
        set jobs($id,callback) $callback
        set jobs($id,stdout) ""
        set jobs($id,tag) $tag

        fileevent $pipe readable [list workerpool::on_readable $id]
    }

    # --------------------------
    # IO Handler
    # --------------------------
    proc on_readable {id} {
        variable jobs

        if {![info exists jobs($id,pipe)]} {
            return
        }

        set pipe $jobs($id,pipe)

        append jobs($id,stdout) [read $pipe]

        if {[eof $pipe]} {
            finalize_job $id
        }
    }

    # --------------------------
    # Finalize Job (Correct Exit Handling)
    # --------------------------
    proc finalize_job {id} {
        variable jobs
        variable running
        variable idle_callback

        set pipe $jobs($id,pipe)

        fileevent $pipe readable {}

        catch {append jobs($id,stdout) [read $pipe]}

        set status [catch {close $pipe} result options]

        if {$status == 0} {
            set exitcode 0
        } else {
            if {[dict exists $options -errorcode]} {
                set ec [dict get $options -errorcode]
                if {[lindex $ec 0] eq "CHILDSTATUS"} {
                    set exitcode [lindex $ec 2]
                } else {
                    set exitcode -1
                }
            } else {
                set exitcode -1
            }
        }

        set cmd $jobs($id,cmd)
        set callback $jobs($id,callback)
        set stdout $jobs($id,stdout)
        set tag $jobs($id,tag)

        foreach key {pipe cmd callback stdout tag} {
            unset jobs($id,$key)
        }

        incr running -1

        # Callback now includes tag
        uplevel #0 [list $callback $id $cmd $exitcode $stdout $tag]

        if {$idle_callback ne ""} {
            uplevel #0 $idle_callback
        }

        check_done
    }

    # --------------------------
    # Cancel Single Job
    # --------------------------
    proc cancel {id} {
        variable queue
        variable jobs

        set newqueue {}
        foreach job $queue {
            if {[lindex $job 0] != $id} {
                lappend newqueue $job
            }
        }
        set queue $newqueue

        if {[info exists jobs($id,pipe)]} {
            foreach p [pid $jobs($id,pipe)] {
                if {$::tcl_platform(platform) eq "windows"} {
                    catch {exec taskkill /PID $p /F}
                } else {
                    catch {exec kill -TERM $p}
                }
            }
        }
    }

    proc cancel_all {} {
        variable queue
        variable jobs

        set queue {}

        foreach key [array names jobs *,pipe] {
            set id [lindex [split $key ,] 0]
            cancel $id
        }

        check_done
    }

    # --------------------------
    # Shutdown Detection
    # --------------------------
    proc check_done {} {
        variable running
        variable queue
        variable done

        if {$running == 0 && [llength $queue] == 0} {
            set done 1
        }
    }

    proc wait {} {
        vwait ::workerpool::done
    }
}



# ---------------------------------------------
# Index Streamer
# ---------------------------------------------

# A data structure to write the energy values in order of frame index
# Since energy values can come in any order due to asyncronous namd calls in order
namespace eval IndexStreamBuffer {

    # Configuration
    variable requiredCount
    variable outputFile

    # Internal state
    variable data
    variable nextIndex 0
    variable fh ""

    array set data {}

    # ------------------------------------------------------------
    # configure <count> <filename>
    # ------------------------------------------------------------
    proc configure {count filename} {
        variable requiredCount
        variable outputFile
        variable fh
        variable data
        variable nextIndex

        if {![string is integer -strict $count] || $count <= 0} {
            error "Count must be a positive integer"
        }

        if {$filename eq ""} {
            error "Filename cannot be empty"
        }

        # Reset state
        array unset data
        set nextIndex 0

        set requiredCount $count
        set outputFile $filename

        # Close old file if open
        if {$fh ne ""} {
            close $fh
            set fh ""
        }

        # Delete existing file
        if {[file exists $outputFile]} {
            file delete -force $outputFile
        }

        # Open new file
        set fh [open $outputFile "w"]
        fconfigure $fh -buffering full
    }

    # ------------------------------------------------------------
    # insert <index> <value>
    # index >= 0
    # ------------------------------------------------------------
    proc insert {index value} {
        variable data
        variable nextIndex

        if {![string is integer -strict $index] || $index < 0} {
            error "Index must be a non-negative integer: given \"$index\""
        }

        set data($index) $value

        flushReady
    }

    proc insert_lines_with_prefix {index src_string prefixes {output_line_suffix ""}} {
		set appendedCount 0;
		set content $src_string;

		set result_str "";

		while {[string length $content] > 0} {
			set idx [string first "\n" $content]

			if {$idx == -1} {
				set line $content
				set content ""
			} else {
				set line [string range $content 0 [expr {$idx - 1}]]
				set content [string range $content [expr {$idx + 1}] end]
			}

			set trimmedLine [string trimleft $line]

			foreach prefix $prefixes {
				set plen [string length $prefix]

				if {[string length $trimmedLine] >= $plen &&
					[string range $trimmedLine 0 [expr {$plen - 1}]] eq $prefix} {

					set result_str "${result_str}${line}${output_line_suffix}";

					incr appendedCount
					break
				}
			}
		}

		if { $result_str ne "" } {
			insert $index $result_str;
		}

		return $appendedCount;
	}

    # ------------------------------------------------------------
    # Flush contiguous blocks forward from nextIndex
    # ------------------------------------------------------------
    proc flushReady {} {
        variable data
        variable nextIndex
        variable requiredCount
        variable fh

        while {1} {

            set endIndex [expr {$nextIndex + $requiredCount}]

            # Check if full block exists
            for {set i $nextIndex} {$i <= $endIndex} {incr i} {
                if {![info exists data($i)]} {
                    return
                }
            }

            # Write block
            for {set i $nextIndex} {$i <= $endIndex} {incr i} {
                puts $fh $data($i)
                unset data($i)
            }

            # Advance window
            set nextIndex [expr {$endIndex + 1}]
        }
    }

    # ------------------------------------------------------------
    # close
    # Flush everything remaining in ascending order
    # ------------------------------------------------------------
    proc close_file {} {
        variable data
        variable nextIndex
        variable fh

        if {$fh eq ""} {
            return
        }

        # Write all remaining indices >= nextIndex in sorted order
        if {[array size data] > 0} {
            set remaining [lsort -integer [array names data]]
            foreach idx $remaining {
                if {$idx >= $nextIndex} {
                    puts $fh $data($idx)
                }
            }
        }

        flush $fh
        close $fh
        set fh ""

        # Clean memory
        array unset data
    }
}


# -------------------------------------------------------------------------------------------



# Helper Function
proc get_dcd_frame_count {filename} {

    if {![file exists $filename]} {
        error "File '$filename' does not exist."
    }

    set fid [open $filename r]
    fconfigure $fid -translation binary -encoding binary

    # Read first 16 bytes (enough to detect format)
    set header16 [read $fid 16]

    if {[string length $header16] < 12} {
        close $fid
        error "File too small to be a valid DCD file."
    }

    # Extract first 8 bytes for record marker detection
    set first8 [string range $header16 0 7]

    # --------------------------------------------------------------
    # Detect:
    #   - 32-bit vs 64-bit record marker
    #   - Endianness
    #
    # Standard DCD first record size = 84 (32-bit marker)
    # Extended CHARMM DCD may use 64-bit record markers
    # --------------------------------------------------------------

    # Try 32-bit little endian
    binary scan $first8 i rec32_le

    # Try 32-bit big endian
    binary scan $first8 I rec32_be

    # Try 64-bit little endian
    binary scan $first8 w rec64_le

    # Try 64-bit big endian
    binary scan $first8 W rec64_be

    set marker_size ""
    set scan_format ""
    set offset 0

    if {$rec32_le == 84} {
        set marker_size 4
        set scan_format "i"
        set offset 4
    } elseif {$rec32_be == 84} {
        set marker_size 4
        set scan_format "I"
        set offset 4
    } elseif {$rec64_le == 84} {
        set marker_size 8
        set scan_format "w"
        set offset 8
    } elseif {$rec64_be == 84} {
        set marker_size 8
        set scan_format "W"
        set offset 8
    } else {
        close $fid
        error "Not a recognized DCD file (invalid record marker)."
    }

    # --------------------------------------------------------------
    # Validate CORD signature
    # --------------------------------------------------------------

    seek $fid $offset start
    set signature [read $fid 4]

    if {$signature ne "CORD"} {
        close $fid
        error "Invalid DCD file (missing CORD signature)."
    }

    # --------------------------------------------------------------
    # Read NSET (number of frames)
    # --------------------------------------------------------------

    set nset_raw [read $fid 4]
    binary scan $nset_raw $scan_format nframes

    close $fid

    return $nframes
}








# -------------------------
# File paths
# -------------------------
set tmp_dir 					"/tmp/namd_energy.dynamic";		# or ".". folder to store all temp files
set namd_temp_files_prefix 		[file join $tmp_dir "${out_file_prefix}.namd-temp"];

set namd_final_log_filename 	"${namd_temp_files_prefix}.log";	# gets deleted in cleanup
set out_erg_filename 			"${out_file_prefix}.energy.csv";

# ensure tmp dir
file mkdir $tmp_dir;


# returns full path of temp_file prefix for a specific frame
proc get_frame_temp_file_prefix { frame_index } {
	global namd_temp_files_prefix

	return "${namd_temp_files_prefix}.frame-${frame_index}";
}


#-----------------------------------
#  Writing pdb, dcd and namd_config file for a single frame in trajectory
#-----------------------------------
proc write_interaction_pdb_and_dcd { frame_index } {
	global mol_id all sel1 sel2 has_sel2 selection1 selection2 atom_count_sel2;
	global namd_temp_files_prefix;

	set file_prefix [get_frame_temp_file_prefix $frame_index];

	animate write dcd "${file_prefix}.dcd" waitfor all beg 0 end 0 $mol_id;

	$all frame 0;
# 	set all [atomselect $mol_id "all" frame 0];
	set sel1 [atomselect $mol_id "$selection1" frame 0];
	if { $has_sel2 == 1 } {
		set sel2 [atomselect $mol_id "$selection2" frame 0];
		set atom_count_sel2 [$sel2 num];
	} else {
		set sel2 "";
		set atom_count_sel2 0;
	}

	puts "------------------------------------------------------------"
	puts " => ATOM COUNT: Selection-1 ([$sel1 num]) | Selection-2 ($atom_count_sel2)"
	puts "------------------------------------------------------------"

	# to avoid large coordinates PDB format cannot handle
	#$all moveby [vecinvert [measure center $all]]
	$all set x 0;
	$all set y 0;
	$all set z 0;

	$all set beta 0;
	$all set occupancy 0;

	$sel1 set beta 1;
	if { $has_sel2 == 1 } {
		$sel2 set beta 2;
	}

	$all writepdb "${file_prefix}.pdb";

	# cleanup
	$sel1 delete;
	if { $has_sel2 == 1 } {
		$sel2 delete;
	}
}

# Write for each frame
proc write_namd_conf { frame_index first_time_step } {
	global psf_file param_files has_ext_sys initial_ext_sys
	global cutoff switchdist dielectric temperature has_sel2 has_pme pme_grid_spacing
	global has_drude

	set file_prefix [get_frame_temp_file_prefix $frame_index];

	set namdconf [open "${file_prefix}.namd" w]

	puts $namdconf "#####################################################################"
	puts $namdconf "# NAMD config to calculate energies from trajectory .dcd files"
	puts $namdconf "# (Auto-generated by NAMD Energy plugin by RC)"
	puts $namdconf "#####################################################################"

	puts $namdconf "\n# Structure ----------------------"
	puts $namdconf "structure \t\t [file normalize ${psf_file}];"
	puts $namdconf "coordinates \t\t [file normalize ${file_prefix}.pdb];"
	puts $namdconf "\noutputname \t\t [file normalize ${file_prefix}];"

	puts $namdconf "\n# Force Field ----------------------"
	puts $namdconf "paraTypeCharmm \t\t on;"
	foreach par $param_files {
		puts $namdconf "parameters \t\t [file normalize ${par}];"
	}

	puts $namdconf "\nexclude \t\t scaled1-4;"
	puts $namdconf "cutoff \t\t\t $cutoff;"
	if {$switchdist == "" || $switchdist < 0.0} {
		puts $namdconf "switching off;"
	} else {
		puts $namdconf "switchdist \t\t $switchdist;"
	}
	puts $namdconf "dielectric \t\t $dielectric;"
	puts $namdconf "temperature \t\t ${temperature};"
	if { $temperature == 0} {
		puts $namdconf "COMmotion \t\t yes;"
	}

	if { $has_ext_sys == 1 } {
		puts $namdconf "\n# Periodic Boundary Condition ----------------------"
		puts $namdconf "extendedSystem \t\t [file normalize ${initial_ext_sys}];"

		if {$has_pme == 1} {
		puts $namdconf "PME \t\t\t on;"
		puts $namdconf "PMEGridSpacing \t\t $pme_grid_spacing; \t # PME grid spacing (in A)"
		}
	}

	# Drude force field
	if { $has_drude == 1 } {
		puts $namdconf "\n# Drude Force Field -------------------------"
		puts $namdconf "drude on"
		puts $namdconf "drudeTemp 1"
		puts $namdconf "drudeDamping 20.0"
		puts $namdconf "drudeBondLen 0.25"
		puts $namdconf "drudeHardWall on"
	#    puts $namdconf "LJcorrection yes"
		puts $namdconf "langevin on"
		puts $namdconf "langevinDamping 1"
		puts $namdconf "langevinTemp $temperature"
		puts $namdconf "langevinHydrogen off"
	}

	# Pair interaction
	puts $namdconf "\n# Pair Interaction ------------------------------"
	puts $namdconf "pairInteraction \t\t on;"
	puts $namdconf "pairInteractionFile \t\t [file normalize ${file_prefix}.pdb];"
	puts $namdconf "pairInteractionCol \t\t B; \t # Beta"
	puts $namdconf "pairInteractionGroup1 \t\t 1;"

	if { $has_sel2 == 1 } {
		puts $namdconf "pairInteractionGroup2 \t\t 2;"
	} else {
		puts $namdconf "pairInteractionSelf \t\t on;"
	}

	puts $namdconf "\n# Calculation Energies from the single frame dcd file ------------------------------"
	puts $namdconf "numsteps \t 1;"
	puts $namdconf "set ts $first_time_step; \t # Timestep counter\n";		# Timestep Counter

	puts $namdconf "coorfile open dcd [file normalize ${file_prefix}.dcd];"
	puts $namdconf "while \{ \[coorfile read\] \!= -1 \} \{"
	puts $namdconf "   firstTimestep \$ts;\n   run 0;\n; "
	puts $namdconf "\}\ncoorfile close;\n"

	flush $namdconf
	close $namdconf
}






# =================================
# Running simulation
# =================================

set _time_namd_str [clock format [clock seconds] -format "%Y-%m-%d %H:%M:%S"]

puts "\n ${_time_namd_str}"
puts "================================================================================================="
puts " -> Running NAMD-Energy: \"$label_\""
puts " -> Selection 1: \"$selection1\" | Selection 2: \"$selection2\""
puts " -> Energies: \"$out_energies\""
if { $has_sel2 == 1 } { puts " -> Forces: $keepforce_str" }
puts " -> Output: \"$out_erg_filename\""
puts "=================================================================================================\n"


# clean all frames
animate delete all;
set ts 		$timestep_first;
set ts_inc	[expr $frame_freq * ($frame_skip + 1)];

set stop_job_generation 0;
set error_occured 0;
set error_log_filename "";

set cur_dcd_index -1; 		# current dcd file index
set cur_dcd_frame_count 0;	# total frame count of current dcd file
set idx 0;					# global frame index across all dcd's
set frames_done 0;			# total frames processed, incremented by dcd_frame_count

set frame_log_merge_count 0;	# how many frame log files are merged to final log file
set line_prefixes { "ETITLE:" "ENERGY:" "PAIR INTERACTION:" };	# line prefixes to append to final log



proc init_next_dcd {} {
	global dcd_files cur_dcd_index cur_dcd_frame_count frames_done

	incr cur_dcd_index 1;
	incr frames_done $cur_dcd_frame_count;

	if { $cur_dcd_index < [llength $dcd_files] } {
		# load variables
		set cur_dcd_frame_count [get_dcd_frame_count [lindex $dcd_files $cur_dcd_index]];

		return 1;
	}

	return 0;
}


# Return "" to stop generation
proc generate_next_job {} {
	global dcd_files cur_dcd_index cur_dcd_frame_count
    global idx frames_done namd_cmd_list
    global ts ts_inc timestep_end mol_id

    if { $timestep_end > 0 && $ts > $timestep_end } {
		return "";
	}

    set i [expr $idx - $frames_done];
    if { $cur_dcd_index < 0 || $i == $cur_dcd_frame_count } {
		set ret [init_next_dcd];

		if { $ret == 0 } {
			# All dcd files processed
			return "";
		}

		set i [expr $idx - $frames_done];	# should be 0
    }


    set dcd_file [lindex $dcd_files $cur_dcd_index];

	puts "-----------------------------------------------"
	puts " => FRAME INDEX : $idx (from_this_file: $i)  | TIMESTEP: $ts"
	puts "-----------------------------------------------"

	# Create DCD and interaction-pdb file for this frame
	mol addfile $dcd_file type dcd first $i last $i waitfor all molid $mol_id;

	# write frame files
	write_interaction_pdb_and_dcd $idx;
	animate delete all;			# clean frame from memory

	# write namd config file
	write_namd_conf $idx $ts;
	incr ts $ts_inc;

	# Create namd exec command
	set frame_tmp_file_prefix [get_frame_temp_file_prefix $idx];
	set frame_namd_cmd [linsert $namd_cmd_list end+1 "${frame_tmp_file_prefix}.namd"]

    incr idx;
    return $frame_namd_cmd;
}


proc job_producer {} {
	global idx
    global stop_job_generation

    if { $stop_job_generation == 1 } {
		return;
	}

    # Fill available slots
    while {[workerpool::capacity] > 0} {

		if { $stop_job_generation == 1 } {
			return;
		}

		set cur_idx $idx;
        set job_cmd [generate_next_job];
        if {$job_cmd eq ""} {
			set stop_job_generation 1;
            break
        }

        workerpool::submit $job_cmd on_job_done $cur_idx; # frame index (idx) as tag
    }
}



proc on_job_failed { id cmd exitcode output idx } {
	global workerpool error_occured stop_job_generation
	global namd_temp_files_prefix namd_final_log_filename error_log_filename

	# Stop generating more jobs
	set error_occured 1;
	set stop_job_generation 1;

	# cancel running jobs
	workerpool::cancel_all

	# write error logs to file
	set fidx_temp_prefix [get_frame_temp_file_prefix $idx];
	set fidx_log_file_name "${fidx_temp_prefix}.log";
	set error_log_filename "$fidx_log_file_name";

	set fid [open "$fidx_log_file_name" "w"];
	puts $fid $output;
	close $fid;

	# Print Log details
	puts "\n-----------------------------------------------------------------";
	puts " => ERROR: NAMD Job Failed for frame index: ${idx}"
	puts " => Final Log file: \"${namd_final_log_filename}\""
	puts " See \"${fidx_log_file_name}\" for error details "
	puts "-----------------------"
	puts "# Cleanup TEMP files after evalaution by running..."
	puts " => rm ${namd_temp_files_prefix}*"
	puts "-------------------------------------------------------------------\n";
}


# Error code is unreliable, so checking for valid output
proc on_job_done {id cmd exitcode output idx} {
	global frame_log_merge_count line_prefixes error_occured IndexStreamBuffer

	# ignore if error occured
	if { $error_occured != 0 } {
		return
	}

	set fidx_temp_prefix [get_frame_temp_file_prefix $idx];
	# set fidx_log_file_name "${fidx_temp_prefix}.log";
	set energy_line_found 0;

	if { [string trim $output] ne "" } {

		# append only energy lines
		set append_count [IndexStreamBuffer::insert_lines_with_prefix $idx $output $line_prefixes "\n\n"];

		if { $frame_log_merge_count == 0 } {
			# Remove ETITLE
			set line_prefixes [lsearch -inline -all -not -exact $line_prefixes "ETITLE:"];
			set energy_line_found [expr $append_count > 1];
		} else {
			set energy_line_found [expr $append_count > 0];
		}

		incr frame_log_merge_count;
	}


	# Success
	if { $energy_line_found > 0 } {
		# Cleanup
		file delete {*}[glob -nocomplain "${fidx_temp_prefix}.*"];
	} else {
		# Failure
		on_job_failed $id $cmd $exitcode $output $idx
	}

}




#-------------------------
# Setup Worker Pool
# ------------------------

set worker_pool_size [expr $namd_processes / $namd_pes];

workerpool::init $worker_pool_size 0;	# no shell support required, commands are given as list

# Whenever a job finishes and a slot frees up, producer gets called automatically.
workerpool::on_idle job_producer


#-------------------------
# Setup IndexStreamBuffer
# ------------------------
set index_buffer_size [expr $worker_pool_size * 2];

IndexStreamBuffer::configure $index_buffer_size $namd_final_log_filename;


#-------------------------
# Start Run
#-------------------------
# initial kickstart
job_producer

# Wait for all jobs to finishe
workerpool::wait



#-------------------------
# Closeup
#-------------------------
# Close buffer (final log file)
IndexStreamBuffer::close_file

$all delete;
mol delete $mol_id




# --------------------------------------------------------------------------------
# Parsing Energies from .log file
# --------------------------------------------------------------------------------

proc listsum {mylist} {
    set sum 0;
	foreach element $mylist {
		set sum [expr $sum + $element];
	}

    return $sum;
}

proc cleanlist {mylist} {
	set newlist [list];

    #Returns a copy of mylist with all empty elements removed
    for {set i 0} {$i < [llength $mylist]} {incr i} {
    		if {[lindex $mylist $i] != {}} {
    			lappend newlist [lindex $mylist $i];
    		}
    }

    return $newlist;
}


set cur_frame 0;
set skip1 [expr $frame_skip + 1];
#set stride [expr $skip1 * $frame_freq]

#Read the input
set log_file [open $namd_final_log_filename "r"];
#     ETITLE: TS BOND ANGLE DIHED IMPRP ELECT VDW BOUNDARY MISC KINETIC TOTAL TEMP POTENTIAL ...
# Index  0    1   2     3     4     5     6    7    8       9     10       11   12     13	    ..

# Output energy file
set fout [open $out_erg_filename "w"];

# Comments ----------------
if { [info exists comment_token] == 1 && [string trim $comment_token] ne "" } {
	set comments [list];

	lappend comments "================ $label_ ================"
	lappend comments "PARAM File(s): \[$param_files\]  (Drude: ${drude})"
	lappend comments "PSF File     : $psf_file"
	lappend comments "DCD File(s)  : \[$dcd_files\]"
	lappend comments "## Selections ---------------"
	lappend comments "SELECTION 1  : \"$selection1\""
	lappend comments "SELECTION 2  : \"$selection2\""
	lappend comments "INITIAL ATOM COUNT   : Total ($atom_count_total) | Selection 1 ($atom_count_sel1) | Selection 2 ($atom_count_sel2)"
	lappend comments "## Output -------------------"
	lappend comments "OUTPUT Energies : $out_energies"
	lappend comments "OUTPUT Force    : $keepforce_str  (Force on Sel-1 due to Sel-2)"
	lappend comments "TIMESTEP First  : $timestep_first \t\t (only for book-keeping)"
	lappend comments "FRAME Frequency : $frame_freq steps \t (only for book-keeping)"
	lappend comments "FRAME Skip      : $frame_skip frames"
	lappend comments "## Simulation Params---------"
	lappend comments "Temperature : $temperature K (KBT = [format %.2e [expr 0.0019872 * $temperature]] kcal/mol = [format %.2e [expr 1.380649e-23 * $temperature]] J/molecule)"
	lappend comments "Cutoff      : $cutoff Å"
	lappend comments "Switchdist  : $switchdist Å"
	lappend comments "Dielectric  : $dielectric"
	lappend comments "## Periodic BC --------------"
	lappend comments "PERIODIC BC    : $has_ext_sys"
	lappend comments "PME            : $pme"
	lappend comments "PME Grid Space : $pme_grid_spacing Å"
	lappend comments "Initial Extended System File : \"$initial_ext_sys\""
	lappend comments "-----------------------------------------------------------------"
	lappend comments "Units => ENERGY: 1 kcal/mol     = 6.95e-21 J/molecule = [format %.2e [expr 1/(0.0019872 * $temperature)]] KBT"
	lappend comments "      => FORCE : 1 kcal/(mol Å) = 69.5 pN"
	lappend comments "========================================================================="

	foreach c $comments {
		puts $fout "${comment_token} ${c}";
	}
}


#Initialize variables to hold lists of the output from all frames
set outputlists [list];

#Set up a header line for the output
#Create the output headings and vector
lappend lhead "FRAME";
lappend lhead "TS";
if ([regexp {all|bond|conf} $out_energies]) {lappend lhead "BOND"};
if ([regexp {all|angl|conf} $out_energies]) {lappend lhead "ANGLE"};
if ([regexp {all|dihe|conf} $out_energies]) {lappend lhead "DHIED"};
if ([regexp {all|impr|conf} $out_energies]) {lappend lhead "IMPRP"};
if ([regexp {all|conf} $out_energies]) {lappend lhead "CONF"};
if ([regexp {all|elec|nonb} $out_energies]) {lappend lhead "ELECT"};
if ([regexp {all|vdw|nonb} $out_energies]) {lappend lhead "VDW"};
if ([regexp {all|nonb} $out_energies]) {lappend lhead "NONBOND"};
if ([regexp {all|boun} $out_energies]) {lappend lhead "BOUNDARY"};
if ([regexp {all|misc} $out_energies]) {lappend lhead "MISC"};
if ([regexp {all|kine} $out_energies]) {lappend lhead "KINETIC"};
if ([regexp {all|pote} $out_energies]) {lappend lhead "POTENTIAL"};
lappend lhead "TOTAL";
if { $keepforce == 1 && [regexp {all|vdw|nonb} $out_energies]} { lappend lhead "VDW_FORCE" };
if { $keepforce == 1 && [regexp {all|elec|nonb} $out_energies]} { lappend lhead "ELECT_FORCE" };
if { $keepforce == 1 } { lappend lhead "TOTAL_FORCE" };

# Header
set headerstring [join $lhead $out_delimiter];
puts $fout "$headerstring";

while {[gets $log_file enerstring] >= 0} {
    #Skip unless we're on an energy line
    if {[regexp {^ENERGY:} $enerstring] == 0} {continue};

    #Initialize the variables to hold data from this frame
    set outputlist [list];
    #Next, make a list with all of the energy fields
    set enerlist [split $enerstring];
    set enerlist [cleanlist $enerlist];

#     ETITLE: TS BOND ANGLE DIHED IMPRP ELECT VDW BOUNDARY MISC KINETIC TOTAL TEMP POTENTIAL ...
# Index  0    1   2     3     4     5     6    7    8       9     10       11   12     13	    ...

    lappend outputlist $cur_frame;			 # Frame
    lappend outputlist [lindex $enerlist 1];	 # Time Step

    ### Confirmational Energies ----------------
    if ([regexp {all|bond|conf} $out_energies]) {lappend outputlist [lindex $enerlist 2] };
    if ([regexp {all|angl|conf} $out_energies]) {lappend outputlist [lindex $enerlist 3] };
    if ([regexp {all|dihe|conf} $out_energies]) {lappend outputlist [lindex $enerlist 4] };
    if ([regexp {all|impr|conf} $out_energies]) {lappend outputlist [lindex $enerlist 5] };

    set conf_energy [listsum [lrange $enerlist 2 5]];  # confimational energy (-conf) = bond + angle + dihed + imprp
    if ([regexp {all|conf} $out_energies]) {lappend outputlist $conf_energy};

    ### Non Bonded Energies ------------------
    if ([regexp {all|elec|nonb} $out_energies]) {lappend outputlist [lindex $enerlist 6] };
    if ([regexp {all|vdw|nonb} $out_energies]) {lappend outputlist [lindex $enerlist 7] };

    set nonb_energy [listsum [lrange $enerlist 6 7]];  # non-bonded energy (-nonb) = elect + vdw
    if ([regexp {all|nonb} $out_energies]) {lappend outputlist $nonb_energy};

    # Other Energies -------------------
    if ([regexp {all|boun} $out_energies]) {lappend outputlist [lindex $enerlist 8] };	# Boundary Energy
    if ([regexp {all|misc} $out_energies]) {lappend outputlist [lindex $enerlist 9] };	# Misc Energy

    # Kinetic Energy
    if ([regexp {all|kine} $out_energies]) {lappend outputlist [lindex $enerlist 10] };

    # Potential Energy = bond + angle + dihed + imprp + elect + vdw + boundary + misc
    if ([regexp {all|pote} $out_energies]) {lappend outputlist [lindex $enerlist 13] };

    # Total Energy = kinetic + potential
    set total_energy [lindex $enerlist 11];				# total = kinetic + potential (-conf + -nonb)
    lappend outputlist $total_energy;

    # Force between selections
    if {$keepforce == 1} {
        gets $log_file;
        gets $log_file forcestring;

        set forcelist [split $forcestring];
        set forcelist [cleanlist $forcelist];

        # Vdw Force ----------
        set forceindex [lsearch -regexp  $forcelist "^VDW_FORCE"];
		set vdwvec [list];
		lappend vdwvec [lindex $forcelist [expr $forceindex + 1]];
		lappend vdwvec [lindex $forcelist [expr $forceindex + 2]];
		lappend vdwvec [lindex $forcelist [expr $forceindex + 3]];

		set vdwproj 1;
		set vdwmag [expr $vdwproj * [veclength $vdwvec]];

		# Electrostatic Force ---------------
		set forceindex [lsearch -regexp  $forcelist "^ELECT_FORCE"];
		set elvec [list];
		lappend elvec [lindex $forcelist [expr $forceindex + 1]];
		lappend elvec [lindex $forcelist [expr $forceindex + 2]];
		lappend elvec [lindex $forcelist [expr $forceindex + 3]];
		set elecproj 1;
		set elecmag [expr $elecproj * [veclength $elvec]];

		# Total Force
		set totforce [expr $vdwmag + $elecmag];

		# Add to output
        if {[regexp {all|vdw|nonb} $out_energies]} {
			lappend outputlist $vdwmag;
        }

        if {[regexp {all|elec|nonb} $out_energies]} {
        		lappend outputlist $elecmag;
        	}

        lappend outputlist $totforce;
    }

	# Formatting Output
    set out_str "[lindex $outputlist 0]${out_delimiter}[lindex $outputlist 1]";		# Frame, TimeStep

    for {set i 2} { $i < [llength $outputlist]} {incr i} {
    		set val_str [format $out_energy_format [lindex $outputlist $i]];
    		set out_str "${out_str}${out_delimiter}${val_str}";
    }

    # Write entry to file
    puts $fout $out_str;

    incr cur_frame $skip1;
    #incr cur_ts $stride;
}


close $log_file;
close $fout;



# ------------------------------------------------------------------------
proc cleanup { } {
	global namd_temp_files_prefix namd_final_log_filename;

    file delete {*}[glob -nocomplain ${namd_temp_files_prefix}*];
    file delete $namd_final_log_filename;
}


# Celaning Up
if { $debug == 0 && $error_occured == 0 } {
	puts "\n--------------------------------------"
	puts " => Deleting temp files: ${namd_temp_files_prefix}*";
	puts "--------------------------------------\n"
	cleanup;
}



set time_end [clock seconds];
set time_end_str [clock format $time_end -format "%Y-%m-%d %H:%M:%S"]

puts "\n ${time_end_str}"
puts "==================  NAMD Energy FINISHED  ====================="
puts "=> LABEL: \"${label_}\""
puts "=> OUTPUT Energy File: \"${out_erg_filename}\""
puts "-> Time Taken: [expr $time_end -$time_start] secs"
if { $error_occured != 0 } {
	puts "-----------------------------------"
	puts "--------- NAMD RUN FAILED ---------"
	puts "-----------------------------------"
	puts " => Final Log file: \"${namd_final_log_filename}\""
	puts " See \"${error_log_filename}\" for error details "
	puts "-----------------------"
	puts "# Cleanup TEMP files after evalaution by running..."
	puts " => rm ${namd_temp_files_prefix}*"
}
puts "===================================================\n"

exit;
