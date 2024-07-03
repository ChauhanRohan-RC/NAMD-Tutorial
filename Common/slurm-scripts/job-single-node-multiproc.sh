#!/bin/bash
#===================================================================
# SLURM Script to run NAMD2 on SINGLE NODE with multiple-proceses
#===================================================================
# TODO: 1. set --job-name, --ntasks-per-node
#	2. set WORK_DIR
#	3. set NAMD .conf and .log files
#	4. submit with "sbatch job-script.sh"  or  "sbatch -w cn03 job-script.sh" 
#
# Job name
#SBATCH --job-name=dna-test-single-node-multiproc
#
# Account
##SBATCH --account=parbati
#
# Partition
##SBATCH --partition=partition_name
#
# Number of nodes
#SBATCH --nodes=1
#
# Tasks per node (based on number of cores per node = 48 in this HPC)
#SBATCH --ntasks-per-node=46
#
# Processors (cores) per task 
#SBATCH --cpus-per-task=1
#
# Standard Output and Error
#SBATCH --output=output.log
#SBATCH --error=error.log
#
# Wall clock Limit
##SBATCH --time=00:00:30
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
WORK_DIR=/scratch/home/parbati/rohan/Projects/MD/dna-hairpin/test
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
module use /scratch/apps/modules
module load openmpi/4.1.6

### ---------------- MAIN SCRIPT ------------------
## APP path
NAMD_HOME=/scratch/home/parbati/NAMD/NAMD_2.14_Linux-x86_64-verbs/ 

## ======= MAIN COMMAND ========

TOTAL_TASKS=$(($SLURM_JOB_NUM_NODES * $SLURM_NTASKS_PER_NODE))

## -> can use +isomalloc_sync option (after namd2)
$NAMD_HOME/charmrun +p$TOTAL_TASKS ++mpiexec ++nodelist $nodefile $NAMD_HOME/namd2 $WORK_DIR/test.conf > $WORK_DIR/test2.log
