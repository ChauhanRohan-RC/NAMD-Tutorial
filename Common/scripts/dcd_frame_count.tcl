#!/bin/tclsh

# Simple Script to get frame count of a DCD file

# Usage:
# ./dcd_frame_count "<filename>.dcd"
#



# ------------------------------------------------------------------
# Function: get_dcd_frame_count
#
# Purpose:
#   Read number of frames (NSET) from a DCD file.
#   Supports:
#       - Standard DCD
#       - CHARMM extended DCD
#       - 32-bit record markers
#       - 64-bit record markers
#       - Auto endianness detection
#
# Arguments:
#   filename - Path to DCD file
#
# Returns:
#   Integer frame count
# ------------------------------------------------------------------

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



# ------------------------------------------------------------------
# Main program
# ------------------------------------------------------------------

# Check command line argument
if {$argc != 1} {
    puts "Usage: tclsh get_dcd_frames.tcl trajectory.dcd"
    exit 1
}

# Get filename from argument
set dcdfile [lindex $argv 0]

# Verify file exists
if {![file exists $dcdfile]} {
    puts "Error: File '$dcdfile' not found."
    exit 1
}

# Get frame count
set frame_count [get_dcd_frame_count $dcdfile]

# Print result
# puts "DCD file: $dcdfile"
# puts "Number of frames: $frame_count"
puts "$frame_count"
