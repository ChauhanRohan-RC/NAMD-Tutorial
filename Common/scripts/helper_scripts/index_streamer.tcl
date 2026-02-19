#!/usr/bin/env tclsh


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
            error "Index must be a non-negative integer"
        }

        set data($index) $value

        flushReady
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



IndexStreamBuffer::configure 4 "test-output.txt"

IndexStreamBuffer::insert 1 "hello 1\n\nsecondline\n\n"
IndexStreamBuffer::insert 4 "hello 4\n\nsecondline\n\n"
IndexStreamBuffer::insert 3 "hello 3\n\nsecondline\n\n"
IndexStreamBuffer::insert 5 "hello 5\n\nsecondline\n\n"
IndexStreamBuffer::insert 2 "hello 2\n\nsecondline\n\n"
IndexStreamBuffer::insert 7 "hello 7\n\nsecondline\n\n"
IndexStreamBuffer::insert 0 "hello 0\n\nsecondline\n\n"

IndexStreamBuffer::close_file
