#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc_%A.out
#SBATCH --error=fastqc_%A.err
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --partition=q_cn

# Description:  QC of raw reads by FASTQC

module load fastqc/0.11.7

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
infile=$1

# Copy necessary files to scratch directory
cp $infile ${SCRATCH_DIRECTORY}
infile_basename=`basename $infile`

# QC by FASTQC
mkdir fastqc
fastqc -t 10 -o fastqc $infile_basename

# After the job is done we copy our output back to our directory
outdir=/home/Magda_lab/maitenat/scratch60/ont
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/fastqc/* ${outdir}