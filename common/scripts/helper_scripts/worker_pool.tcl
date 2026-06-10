#!/usr/bin/env tclsh

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







# -----------------------------
# MAIN
# -----------------------------

set job_counter 0
# set max_generate 20     ;# simulate unknown upper bound
# set stop_generation 0

# Return "" to stop generation
proc generate_job {} {
    global job_counter
#     global max_generate stop_generation

#     if {$job_counter >= $max_generate} {
#         set stop_generation 1
#         return ""
#     }

    if { $job_counter >= 10 } {
        return "";  # done
    }

    incr job_counter
#     return [list namd3 test.conf];    # works without sheel, efficient
    return "sleep 2 && echo hello && exit 3";  # must enable shell
}

proc job_done {id_ cmd exitcode output tag_} {
    puts "FINISHED => cmd: $cmd | Exit code: $exitcode | Output: $output | tag: ${tag_}"
}

proc producer {} {
#     global stop_generation

    # Fill available slots
    while {[workerpool::capacity] > 0} {

#         if { $stop_generation == 1 } {
#             break
#         }
#
        set cmd [generate_job]
        if {$cmd eq ""} {
            break
        }
        workerpool::submit $cmd job_done "TAGME-[clock seconds]"
    }
}


# workerpool::init 3 0;     # no-shell support, efficient
workerpool::init 4 1;       # enable shell support, for shell operators like ">" or "&&"

# Whenever a job finishes and a slot frees up,
# producer gets called automatically.
workerpool::on_idle producer

# Initial kickstart
producer

# Wait until all running jobs complete
workerpool::wait

puts "All dynamic jobs completed."

