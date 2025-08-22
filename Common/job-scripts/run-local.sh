#!/bin/bash

# Specify NAMD directory
NAMD_HOME=$(dirname $(which namd3))
TOTAL_TASKS=2

$NAMD_HOME/namd3 +p$TOTAL_TASKS dna_wb_eq.conf > dna_wb_eq.log
echo "DONE 0"
#$NAMD_HOME/namd3 +p$TOTAL_TASKS dna_wb_eq1.conf > dna_wb_eq1.log
#echo "DONE 1"
#$NAMD_HOME/namd3 +p$TOTAL_TASKS dna_gbis_eq2.conf > dna_gbis_eq2.log
#echo "DONE 2"
#$NAMD_HOME/namd3 +p$TOTAL_TASKS dna_gbis_eq3.conf > dna_gbis_eq3.log
#echo "DONE 3"
#$NAMD_HOME/namd3 +p$TOTAL_TASKS dna_gbis_eq4.conf > dna_gbis_eq4.log
#echo "DONE 4"
#$NAMD_HOME/namd3 +p$TOTAL_TASKS dna_gbis_eq5.conf > dna_gbis_eq5.log
#echo "DONE 5"

# ------------------------------------------------------------
# MAIL ON COMPLETION
host="$HOSTNAME"
sender="${host^^} <pbiswaslab@gmail.com>"
receiver="rsingh1.phd@chemistry.du.ac.in"

subject="$host - ${PWD##*/} DONE"
message=" $host : NAMD ${PWD##*/} completed with exit code: $?
\n\n--------------\n# LOGS \n-------------\n\
$(tail -n 15 *.log)"

echo -e "$message" | mail -r "$sender" -s "$subject" "$receiver"
# ------------------------------------------------------------
