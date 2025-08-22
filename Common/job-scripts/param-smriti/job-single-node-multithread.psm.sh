#!/bin/bash
#
#===================================================================
# SLURM Script to run NAMD2 on SINGLE NODE with multiple-Threads in a single process
#===================================================================
# TODO: 1. set --job-name, --ntasks-per-node, --time
#	2. set WORK_DIR
#	3. set NAMD .conf and .log files
#	4. submit with "sbatch job_script.sh"  or  "sbatch -w cn03 job_script.sh" 
#
# Job name
#SBATCH --job-name=dna-t3
#
# Account
##SBATCH --account=parbatib
#
# Partition
#SBATCH --partition=standard
#
# Number of nodes
#SBATCH --nodes=1
#
# Number of gpu's to use --gres=gpu:<num_gpu_to_use>
##SBATCH --gres=gpu:1
#
# Tasks per node (based on number of cores per node = $SLURM_CPUS_ON_NODE)
#SBATCH --ntasks-per-node=48
#
# Processors (cores) per task 
#SBATCH --cpus-per-task=1
#
# Wall clock Time Limit (days-hr:min:secs)
#SBATCH --time=1-00:00:00
#
# Standard Output and Error
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

## Working Dir: Mostly $SLURM_SUBMIT_DIR
WORK_DIR=$SLURM_SUBMIT_DIR
cd $WORK_DIR

### -------- Creating Nodelist file for NAMD ---------
nodes=$(scontrol show hostname $SLURM_JOB_NODELIST)

nodefile=${WORK_DIR}/${SLURM_JOB_NAME}.nodelist

echo "group main" >| $nodefile 
for n in ${nodes[@]}; do
	echo "host ${n}" >> $nodefile
done
### ------------------------------------------------

## Module Load - MPI (NO NEED FOR SINGLE NODE MULTITHREADED MODE)
##module use /scratch/apps/modules
#module load openmpi3/3.1.4

### ---------------- MAIN SCRIPT ------------------
## APP path (Multicore Version)
NAMD_HOME=/home/parbatib/NAMD/NAMD_2.14_Linux-x86_64-multicore

## ======= MAIN COMMAND ========
TOTAL_TASKS=$(($SLURM_JOB_NUM_NODES * $SLURM_NTASKS_PER_NODE))

## -> can use +isomalloc_sync option (after namd2)
$NAMD_HOME/namd2 +p$TOTAL_TASKS $WORK_DIR/dna_gbis_pcf.conf > $WORK_DIR/dna_gbis_pcf.log
