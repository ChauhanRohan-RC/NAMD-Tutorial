#!/bin/bash
#===================================================================
# SLURM Script to run NAMD2 on Multiple Nodes using open-mpi
#===================================================================
# TODO: 1. set --job-name, --nodes, --ntasks-per-node, --time
#	2. set WORK_DIR
#	3. set NAMD .conf and .log files 
#	4. submit with "sbatch job_script.sh"  or  "sbatch -w cn[01-03] job_script.sh" 
#
# Job name
#SBATCH --job-name=dna-t1
#
# Account
##SBATCH --account=parbatib
#
# Partition
#SBATCH --partition=standard
#
# Number of nodes
#SBATCH --nodes=4
#
## Number of gpu's to use --gres=gpu:<num_gpu_to_use>
##SBATCH --gres=gpu:1
#
# Tasks (processes) per node (based on number of cores per node = 48 in this HPC)
#SBATCH --ntasks-per-node=48
#
# Processors (cores) per task 
#SBATCH --cpus-per-task=1
#
# Wall clock Time Limit (days-hr:min:secs)
#SBATCH --time=03:00:00
#
# Standard Output and Error
#SBATCH --output=output.log
#SBATCH --error=error.log
#
##SBATCH --mail-type=ALL
##SBATCH --mail-user=rsingh1.phd@chemistry.du.ac.in

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

## Module Load - MPI
#module use /scratch/apps/modules
module load openmpi3/3.1.4

### ---------------- MAIN SCRIPT ------------------
## APP path
NAMD_HOME=/home/parbatib/NAMD/NAMD_2.14_Linux-x86_64-verbs

## ======= MAIN COMMAND ========

TOTAL_TASKS=$(($SLURM_JOB_NUM_NODES * $SLURM_NTASKS_PER_NODE))

## -> can use +isomalloc_sync option (after namd2)
$NAMD_HOME/charmrun +p$TOTAL_TASKS ++mpiexec ++nodelist $nodefile $NAMD_HOME/namd2 $WORK_DIR/dna_gbis_pcf.conf > $WORK_DIR/dna_gbis_pcf.log
