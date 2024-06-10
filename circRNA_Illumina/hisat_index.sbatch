#!/bin/bash
#SBATCH --job-name=hisatidx
#SBATCH --output=hisat-idx_%A.out
#SBATCH --error=hisat-idx_%A.err
#SBATCH --ntasks=1
#SBATCH --partition=q_fat
#SBATCH --mem-per-cpu=30G

# Description: Index the genome with hisat2

module load hisat2/2.1.0

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}

genome=/home/Magda_lab/maitenat/DATA/indeces/mm39.fa

hisat2-build -p 4 $genome mm39_ht2