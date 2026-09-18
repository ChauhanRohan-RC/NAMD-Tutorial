#!/usr/bin/env tclsh


# A data structure to write the energy values in order of frame index
# Since energy values can come in any order due to asyncronous namd calls
namespace eval IndexStreamBuffer {

    # Configuration
    variable out_file_path;
    variable chunk_size;      # number of enties in a chunk

    # Internal state
    variable data;
    variable next_index 0;
    variable fh "";
    variable written_chunk_count 0;

    array set data {};

    # ------------------------------------------------------------
    # configure <chunk_size> <output_file_path>
    # ------------------------------------------------------------
    proc configure {required_chunk_size output_file_path} {
        variable chunk_size
        variable out_file_path
        variable fh
        variable data
        variable next_index
        variable written_chunk_count

        if {![string is integer -strict $required_chunk_size] || $required_chunk_size <= 0} {
            error "chunk size must be a positive integer"
        }

        if {$output_file_path eq ""} {
            error "output file path cannot be empty"
        }

        # Reset state
        array unset data
        set next_index 0
        set written_chunk_count 0

        set chunk_size $required_chunk_size
        set out_file_path $output_file_path

        # Close old file if open
        if {$fh ne ""} {
            close $fh
            set fh ""
        }

        # Delete existing file
        if {[file exists $out_file_path]} {
            file delete -force $out_file_path
        }

        # Open new file
        set fh [open $out_file_path "w"]
        fconfigure $fh -buffering full
    }

    # ------------------------------------------------------------
    # insert <index> <value>
    # index >= 0
    # ------------------------------------------------------------
    proc insert {idx value} {
        variable data
        variable next_index

        if {![string is integer -strict $idx] || $idx < 0} {
            error "Index must be a non-negative integer"
        }

        set data($idx) $value

        consider_flush
    }

    # ------------------------------
    # Callbacks
    #-------------------------------

    proc pre_chunk_write { chunk_index } {
        puts "PRE_CHUNK_WRITE: $chunk_index"
    }

    proc post_chunk_write { chunk_index } {
        puts "POST_CHUNK_WRITE: $chunk_index"
    }

    # Convert data entry to output string
    proc convert_to_out_string { data_entry } {
        return $data_entry
    }

    # ------------------------------------------------------------
    # Flush contiguous blocks forward from next_index
    # ------------------------------------------------------------
    proc consider_flush {} {
        variable data
        variable next_index
        variable chunk_size
        variable fh
        variable written_chunk_count

        while {[array size data] >= $chunk_size} {

            set endIndex [expr {$next_index + $chunk_size}]

            # Check if full block exists
            for {set i $next_index} {$i < $endIndex} {incr i} {
                if {![info exists data($i)]} {
                    return
                }
            }

            set chunk_index $written_chunk_count;

            # pre chunk write
            pre_chunk_write $chunk_index;

            # Write block
            for {set i $next_index} {$i < $endIndex} {incr i} {
                puts $fh [convert_to_out_string $data($i)];
                unset data($i);
            }

            # Advance window
            set next_index $endIndex;
            incr written_chunk_count;

            # post chunk write
            post_chunk_write $chunk_index;
        }
    }

    # ------------------------------------------------------------
    # close
    # Flush everything remaining in ascending order
    # ------------------------------------------------------------
    proc close_file {} {
        variable data
        variable next_index
        variable fh
        variable written_chunk_count

        if {$fh eq ""} {
            return
        }

        # Write all remaining indices >= next_index in sorted order
        if {[array size data] > 0} {
            set remaining [lsort -integer [array names data]]
            set chunk_index $written_chunk_count;

            # pre chunk write
            pre_chunk_write $chunk_index;

            set c 0;
            foreach idx $remaining {
                if {$idx >= $next_index} {
                    puts $fh [convert_to_out_string $data($idx)];
                    incr c;
                }
            }

            # Advance window
            incr next_index $c;
            incr written_chunk_count;

            # post chunk write
            post_chunk_write $chunk_index;
        }

        flush $fh;
        close $fh;
        set fh "";
        set next_index 0;
        set written_chunk_count 0;

        # Clean memory
        array unset data;
    }
}





# -------------------------------------------
# TEST Implementation
# -------------------------------------------

proc shuffleList {list} {
    set len [llength $list]
    for {set i [expr {$len - 1}]} {$i > 0} {incr i -1} {
        # Pick a random index from 0 to i
        set j [expr {int(rand() * ($i + 1))}]

        # Swap elements at index i and j
        set temp [lindex $list $i]
        lset list $i [lindex $list $j]
        lset list $j $temp
    }
    return $list
}


set indices [list];
for {set i 0} {$i < 100} {incr i} {
    lappend indices $i
}

set indices [shuffleList $indices];       # shuffle indices
puts "SHUFFLED INDICES: $indices \n"

IndexStreamBuffer::configure 13 "test_index_streamer.txt"
foreach i $indices {
    IndexStreamBuffer::insert $i "VALUE LINE : $i"
}

IndexStreamBuffer::close_file
