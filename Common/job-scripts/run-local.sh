#!/bin/bash
# Script to run NAMD locally (SINGLE NODE)

# =========================
# INPUT Params
# =========================

## NAMD ---------------
NAMD_DIR=$NAMD_MULTICORE				# or $(dirname $(which namd3))
TOTAL_TASKS=2						# threads TODO
wait_other_namd_processes=false		# [true/false]

namd_cmd="${NAMD_DIR}/namd3 +p${TOTAL_TASKS} +setcpuaffinity"

## Runs ---------------
prefix="dna_wb_eq"			# TODO
error_prefix="error"

# run indices range to execute
run_index_start=0		# [INCLUSIVE] first run index	TODO
run_index_end=1			# [INCLUSIVE] last run index		TODO
run_index_inc=1			# run index increment (DEFAULT: 1)

## Post Run -----------
send_mail=true			# [true/false] Send mail on finish




# =========================
# Execute
# =========================

# wait for other NAMD processes (if they are running)
if $wait_other_namd_processes; then
	waitpid $(pgrep NAMD)
fi

for (( i=$run_index_start; i<=$run_index_end; i+=$run_index_inc )); do
	if [ $i -eq 0 ]; then
		# First Run (Index 0) => No numerical suffix
		echo "=> Starting RUN ${i}"
		$namd_cmd "${prefix}.conf" > "${prefix}.log" 2> "${error_prefix}.log"
		echo "=> RUN ${i} FINISHED"
	elif [ $i -gt 0 ]; then
		# Continuation Run => Use index as suffix
		echo "=> Starting RUN ${i}"
		$namd_cmd "${prefix}${i}.conf" > "${prefix}${i}.log" 2> "${error_prefix}${i}.log"
		echo "=> RUN ${i} FINISHED"
	else
		echo "INVALID RUN INDEX: ${i}. Must be >= 0. Skipping..."
	fi
done


# =========================
# MAIL on Completion
# =========================
host="$HOSTNAME"
sender="${host^^} <pbiswaslab@gmail.com>"
receiver="rsingh1.phd@chemistry.du.ac.in"

do_send_mail() {
	subject="$host - ${PWD##*/} DONE"
	message=" $host : NAMD ${PWD##*/} completed with exit code: $?
	\n\n--------------\n# LOGS \n-------------\n\
	$(tail -n 15 *.log)"

	echo -e "$message" | mail -r "$sender" -s "$subject" "$receiver"
}

if $send_mail; then
	do_send_mail
fi

