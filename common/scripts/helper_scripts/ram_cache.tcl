#!/bin/tclsh

# ---------------------------------------------------
# A RAM Cache to copy files temporarily to RAM
# Mostly used for loading .dcd files to RAM for faster access
#
# call purge to delete all cached files on exit or in exit handlers
# OS Signal caching (HUP INT TERM) requires Tclx extension
#----------------------------------------------------
namespace eval ramcache {
    variable ramdisk "";
    variable buffer_size 2147483648;     # bytes to leave as buffer in ramdisk, default 2GiB
    variable tracked_files [list];
    variable file_map;
    array set file_map {};

    variable net_loading_time_secs 0;    # total time (secs) spent on loading files to ram disk

    proc init {ramdisk_path {buffer_bytes 2147483648}} {
        variable ramdisk
        variable buffer_size
        variable tracked_files
        variable file_map
        variable net_loading_time_secs

        set ramdisk $ramdisk_path
        set buffer_size $buffer_bytes
        set tracked_files [list]
        array unset file_map *
        set net_loading_time_secs 0;

        # ensure temp dir
        file mkdir $ramdisk;
        puts "ramcache: Initialized on $ramdisk with a [expr {$buffer_size / 1024 / 1024}] MiB safety buffer."
    }

    proc get_net_load_time_secs {} {
        variable net_loading_time_secs
        return $net_loading_time_secs;
    }

    # Helper: Get available space in bytes on Linux
    proc get_free_space {dir} {
        set df_out [exec df -B1 $dir]
        set lines [split $df_out "\n"]
        set data_line [lindex $lines 1]

        set fields [regexp -all -inline {\S+} $data_line]
        return [lindex $fields 3]
    }

    # Query the cached RAM path using the original disk path
    proc get_cached_path {orig_path} {
        variable file_map

        # Iterate through the map to find the matching original path
        foreach cached_path [array names file_map] {
            if {$file_map($cached_path) eq $orig_path} {
                return $cached_path
            }
        }

        # Return empty string if not found
        return ""
    }


    # Main function to cache the DCD file
    proc load_file {filepath} {
        variable ramdisk
        variable buffer_size
        variable tracked_files
        variable file_map
        variable net_loading_time_secs

        if {$ramdisk eq ""} {
            error "ramcache Error: Not initialized. Call ::RamCache::init first."
        }
        if {![file exists $filepath]} {
            error "ramcache Error: File $filepath does not exist."
        }

        set fsize [file size $filepath]
        set free_space [get_free_space $ramdisk]
        set required_space [expr {$fsize + $buffer_size}]

        if {$free_space > $required_space} {
            set filename [file tail $filepath]
            set root [file rootname $filename]
            set ext [file extension $filename]

            # Generate a highly unique suffix using process ID and millisecond clock
            set unique_suffix "[pid]_[clock clicks -milliseconds]"
            set new_filename "${root}_${unique_suffix}${ext}"
            set newpath [file join $ramdisk $new_filename]

            # 1. TRACK BEFORE COPYING
            # This ensures that if a signal kills the script during the copy,
            # the signal handler already has the path and will delete the partial file.
            lappend tracked_files $newpath
            set file_map($newpath) $filepath

            # 2. CATCH COPY ERRORS
            # Wrap the blocking copy operation in a catch block
            if {[catch {

                puts "\n---------------------------------------------"
                puts "ramcache: copying file to RAM ...."
                puts " => $filepath --> $newpath"
                puts "\n---------------------------------------------\n"

                set stime [clock seconds];
                file copy -force $filepath $newpath;
                set etime [clock seconds];

                set ttaken [expr $etime - $stime];
                incr net_loading_time_secs $ttaken;

                puts "-----------------------------------------------"
                puts "ramcache: file copied to RAM disk"
                puts " => $filepath --> $newpath"
                puts "ramcache: time taken: $ttaken secs"
                puts "-----------------------------------------------"
            } copy_err]} {
                # If we get here, the copy failed mid-way (e.g., source drive error)
                puts "ramcache ERROR: Failed to copy file to RAM ($copy_err)."
                puts "ramcache: Cleaning up partial file and falling back to HDD..."

                # Delete the partially copied file immediately
                if {[file exists $newpath]} {
                    file delete -force $newpath
                }

                # Remove it from the tracking list so purge doesn't try to delete it later
                set idx [lsearch -exact $tracked_files $newpath]
                if {$idx >= 0} {
                    set tracked_files [lreplace $tracked_files $idx $idx]
                }
                catch {unset file_map($newpath)}

                # Fallback: return the original HDD filepath so the analysis can still proceed
                return $filepath
            }

            # Copy succeeded flawlessly
            return $newpath
        } else {
            puts "ramcache: Insufficient space in $ramdisk. Reading directly from HDD."
            return $filepath
        }
    }


    # Function to safely delete all cached files
    proc purge {} {
        variable tracked_files
        variable file_map

        puts "ramcache: Purging tracked files from RAM..."

        if {[llength $tracked_files] == 0} { return }

        foreach f $tracked_files {
            if {[file exists $f]} {
                file delete -force $f
                puts "ramcache: Deleted $f"
            }
        }
        set tracked_files [list]
        array unset file_map *
    }

    # Delete a specific file from ramcache using its CACHED RAM path
    proc delete_by_cached_path {cached_path} {
        variable tracked_files
        variable file_map

        set idx [lsearch -exact $tracked_files $cached_path]
        if {$idx >= 0} {
            if {[file exists $cached_path]} {
                file delete -force $cached_path
                puts "ramcache: Deleted cached file $cached_path"
            }
            set tracked_files [lreplace $tracked_files $idx $idx]
            catch {unset file_map($cached_path)}
        } else {
            puts "ramcache Warning: $cached_path is not currently tracked."
        }
    }

    # Delete a specific file from ramcache using its ORIGINAL disk path
    proc delete_by_original_path {orig_path} {
        variable file_map
        set found 0

        foreach cached_path [array names file_map] {
            if {$file_map($cached_path) eq $orig_path} {
                delete_by_cached_path $cached_path
                set found 1
            }
        }

        if {!$found} {
            puts "ramcache Warning: No cached file found for original path: $orig_path"
        }
    }
}




# TCL exit handler
proc exit_handler {args} {
    puts "LOG: ON EXIT"
    ramcache::purge
}

# OS signal handler (HUP TERM INT). Requires Tclx extension
proc os_signal_handler {args} {
    puts "LOG: CAUGHT OS SIGNAL $args"

    ramcache::purge
#     exit 1;
}

# Set up the exit/crash handlers for cleanup
proc setup_exit_handlers {} {
    # Catch standard TCL exits
    trace add execution exit enter exit_handler;

    # Try to trap POSIX signals (requires Tclx)
    if {![catch {package require Tclx}]} {
        signal trap {HUP TERM INT} os_signal_handler;
    } else {
        puts "Warning: Tclx package not found. Native OS signal trapping disabled.";
    }
}



setup_exit_handlers

ramcache::init "/tmp/namd_energy.dynamic"
ramcache::load_file "../amyl_wb_eq2.dcd"
