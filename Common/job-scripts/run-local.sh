#!/bin/bash
# Script to run NAMD locally

# =========================
# NAMD Params
# =========================
#NAMD_HOME=$(dirname $(which namd3))
NAMD_DIR=$NAMD_MULTICORE
TOTAL_TASKS=2

# =========================
# Execute
# =========================
$NAMD_DIR/namd3 +p$TOTAL_TASKS dna_wb_eq.conf > dna_wb_eq.log 2> error.log
echo "DONE 0"
$NAMD_DIR/namd3 +p$TOTAL_TASKS dna_wb_eq1.conf > dna_wb_eq1.log 2> error1.log
echo "DONE 1"
$NAMD_DIR/namd3 +p$TOTAL_TASKS dna_wb_eq2.conf > dna_wb_eq2.log 2> error2.log
echo "DONE 2"
#$NAMD_DIR/namd3 +p$TOTAL_TASKS dna_wb_eq3.conf > dna_wb_eq3.log 2> error3.log
#echo "DONE 3"
#$NAMD_DIR/namd3 +p$TOTAL_TASKS dna_wb_eq4.conf > dna_wb_eq4.log 2> error4.log
#echo "DONE 4"
#$NAMD_DIR/namd3 +p$TOTAL_TASKS dna_wb_eq5.conf > dna_wb_eq5.log 2> error5.log
#echo "DONE 5"

# =========================
# MAIL on Completion
# =========================
SEND_MAIL=true

host="$HOSTNAME"
sender="${host^^} <pbiswaslab@gmail.com>"
receiver="rsingh1.phd@chemistry.du.ac.in"

send_mail() {
	subject="$host - ${PWD##*/} DONE"
	message=" $host : NAMD ${PWD##*/} completed with exit code: $?
	\n\n--------------\n# LOGS \n-------------\n\
	$(tail -n 15 *.log)"

	echo -e "$message" | mail -r "$sender" -s "$subject" "$receiver"
}

if $SEND_MAIL; then
	send_mail
fi

