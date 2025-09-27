#!/bin/bash
# =================================================================
# Script to run NAMD locally (SINGLE NODE)
# =================================================================
# Usage: nohup ./run-local.sh > run-local.log & 

# =========================
# INPUT Params
# =========================

## NAMD ---------------
NAMD_DIR=$NAMD_MULTICORE				# or $(dirname $(which namd3))
TOTAL_TASKS=3						# threads TODO
wait_other_namd_processes=true		# [true/false]

namd_cmd="${NAMD_DIR}/namd3 +p${TOTAL_TASKS} +setcpuaffinity"

## Runs ---------------
prefix="amyld_wb_eq"			# TODO
error_prefix="error"

# run indices range to execute
run_index_start=0		# [INCLUSIVE] first run index	TODO
run_index_end=1			# [INCLUSIVE] last run index		TODO
run_index_inc=1			# run index increment (DEFAULT: 1)

## Post Run -----------
send_mail=true			# [true/false] Send mail on finish
require_exit_code0=true	# [true/false] If a RUN returns non-zero exit code, 
						# do not start subsequent RUNS and quit

purge_old_bak_files=true	 # delete .old and .BAK files after runs





# ====================================
# MAIN
# ====================================
run_logs=""
log() {
	run_logs="${run_logs}\n$1"
	echo -e "$1"
}

# Purges unnecessry files
# => USAGE: purge_files <prefix> <error_prefix>
purge_files() {
	local pre="$1"
	local err_pre="$2"

	# delete empty error.log
	if [[ -n "$err_pre" ]]; then
		find "." -name "${err_pre}.log" -type f -empty -delete
	fi

	if [[ -n "$pre" ]]; then
		find "." -name "${pre}.*.old" -type f -delete	# delete *.old files
		find "." -name "${pre}.*.BAK" -type f -delete	# delete *.BAK files
	fi
}


# -------------------------------
# Wait for other NAMD processes
# -------------------------------
if [ $wait_other_namd_processes == "true" ]; then
	while true; do
		namd_pids="$(pgrep NAMD)"
		if [ -z "$namd_pids" ]; then
			break
		fi

		namd_pids="${namd_pids//$'\n'/ }"	# replace "\n" with " "
		log "$(date): Waiting for NAMD processes: ${namd_pids}"
		tstart=$SECONDS
		eval "waitpid $namd_pids"
		tend=$SECONDS
		log "$(date): Done waiting for $((tend - tstart)) secs"
		sleep 2		# Time buffer for other scripts to start subsequent NAMD processes
	done
fi


# --------------------------
# Execute RUNS
# --------------------------
for (( i=$run_index_start; i<=$run_index_end; i+=$run_index_inc )); do
	# Checks ------------------------------
	if [ $i -lt 0 ]; then
		log "ERROR: INVALID RUN INDEX: ${i}. Must be >= 0. Skipping..."
		continue
	fi

	num_suffix="${i}"
	if [ $i -eq 0 ]; then
		num_suffix=""		# First Run (Index 0) => No numerical suffix
	fi

	# Execute NAMD command -----------------
	log "\n$(date): RUN ${i} starting"
	tstart=$SECONDS

	$namd_cmd "${prefix}${num_suffix}.conf" > "${prefix}${num_suffix}.log" 2> "${error_prefix}${num_suffix}.log"
	exit_code="$?"

	tend=$SECONDS
	log "$(date): RUN ${i} finished => Exit Code: ${exit_code} | Duration: $((tend -tstart)) secs"
	# --------------------------------------


	# Purge old, bak and empty logs on successful run
	if [[ $exit_code -eq 0 && $purge_old_bak_files == "true" ]]; then
		purge_files "${prefix}${num_suffix}" "${error_prefix}${num_suffix}"
	fi

	# Sequential Runs: Exit on Failure
	if [[ $exit_code -ne 0 && $require_exit_code0 == "true" ]]; then
		log " => Non-Zero Exit code encountered in RUN ${i}. Quitting..."
		break
	fi
done




# ============================
# MAIL on Completion
# ============================
do_send_mail() {
	host="$(hostname -f)"
	sender="${host^^} <pbiswaslab@gmail.com>"
	receiver="rsingh1.phd@chemistry.du.ac.in"

	runs_str="$(seq ${run_index_start} ${run_index_incr} ${run_index_end})"
	runs_str="${runs_str//$'\n'/ }"

	namd_logs="=====================\n# NAMD LOGS \n====================="
	for (( i=$run_index_start; i<=$run_index_end; i+=$run_index_inc )); do

		if [ $i -eq 0 ]; then
			log_files="${prefix}.log ${error_prefix}.log"
		else
			log_files="${prefix}${i}.log ${error_prefix}${i}.log"
		fi

		log_entry=$(eval "tail ${log_files} 2>/dev/null")
		if [ -n "$log_entry" ]; then
			namd_logs="${namd_logs}\n\n${log_entry}"
		fi
	done


	#subject="$host - ${PWD##*/} DONE"
	subject="${host} - ${prefix} DONE"
	message=" ${host^^} : NAMD ${prefix} RUNS [${runs_str}] Finished\n\n\
=====================\n# RUN LOGS \n=====================\
${run_logs}\n\n\
${namd_logs}"

	echo -e "$message" | mail -r "$sender" -s "$subject" "$receiver"
	postqueue -f
	log "$(date): Mail sent to $receiver"
}

if $send_mail; then
	do_send_mail
fi

