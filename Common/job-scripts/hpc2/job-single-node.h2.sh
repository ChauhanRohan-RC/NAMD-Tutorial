#!/bin/bash
#
#===================================================================
# SLURM Script to run NAMD3 on SINGLE NODE (multicore)
#===================================================================
# TODO: 1. set --job-name, --ntasks-per-node, --time
#	2. set WORK_DIR
#	3. set NAMD .conf and .log files
#	4. submit with "sbatch job_script.sh"  or  "sbatch -w cn03 job_script.sh" 
#
## Job name
#SBATCH --job-name=dna-t3
#
## Account
##SBATCH --account=parbatib
#
## Partition
#SBATCH --partition=cpu
#
## Number of nodes
#SBATCH --nodes=1
#
## Tasks (processes) per node (best is 1 process per node)
#SBATCH --ntasks-per-node=1
#
## CPUs (cores) per process (based on no. of cores in the node = 48 (cpu), 40 (gpu) $SLURM_CPUS_ON_NODE)
#SBATCH --cpus-per-task=48
#
## Number of GPU's to use --gres=gpu:<num_gpu_to_use>
##SBATCH --gres=gpu:2
#
## Wall clock Time Limit (days-hr:min:secs)
#SBATCH --time=1-00:00:00
#
## Standard Output and Error
#SBATCH --output=output.log
#SBATCH --error=error.log
#
##SBATCH --mail-type=ALL
##SBATCH --mail-user=chauhanrohanrc803@gmail.com

##PBS -N dna-test1
##PBS -l nodes=3:ppn=16  
##PBS -l walltime=10:10:00  
##PBS -o output.log      
##PBS -e error.log      
##PBS -V

load_modules() {
	module use "/scratch/home/parbati/programs/modules"

	module load gcc/9.4.0
	#module load openmpi/4.1.6
	module load namd3
	module load vmd
}

load_modules
nodefile="${SLURM_JOB_NAME}.nodelist"


## NAMD ---------------
NAMD_DIR=$NAMD_MULTICORE
TOTAL_TASKS=$(($SLURM_JOB_NUM_NODES * $SLURM_NTASKS_PER_NODE * $SLURM_CPUS_PER_TASK))		# threads TODO

namd_cmd="${NAMD_DIR}/namd3 +p${TOTAL_TASKS} +setcpuaffinity"

## Runs ---------------
work_dir=$SLURM_SUBMIT_DIR		# work dir, mostly $SLURM_SUBMIT_DIR
prefix="amyld_wb_eq"			# TODO
error_prefix="error"

# run indices range to execute
run_index_start=0		# [INCLUSIVE] first run index		TODO
run_index_end=1			# [INCLUSIVE] last run index		TODO
run_index_inc=1			# run index increment (DEFAULT: 1)

## Post Run -----------
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


create_nodelist() {
	nodes=$(scontrol show hostname $SLURM_JOB_NODELIST)
	local filename="$1"

	echo "group main" >| $filename
	for n in ${nodes[@]}; do
		echo "host ${n}" >> $filename
	done
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



# --------------------------
# Execute RUNS
# --------------------------

cd $work_dir
create_nodelist "$nodefile"


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

