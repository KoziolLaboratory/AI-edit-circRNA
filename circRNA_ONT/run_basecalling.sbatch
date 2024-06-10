#!/bin/bash
#SBATCH --job-name=prepNano
#SBATCH --output=preprocess_nanopore_%A.out
#SBATCH --error=preprocess_nanopore_%A.err
#SBATCH -n 1
#SBATCH --partition=q_ai 
#SBATCH --mem-per-cpu=100G
#SBATCH --gres=gpu:1

# Description: Basecalling of Nanopore reads using Guppy

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
indir=$1
indir_basename0=`dirname $indir`
indir_basename=`basename ${indir_basename0}`
mkdir basecalled

module load guppy/6.1.2_gpu

## GPU
guppy_basecaller -i $indir \
-s basecalled \
-x 'cuda:0' \
--config /usr/nzx-cluster/apps/guppy/6.1.2/gpu/data/dna_r10.4_e8.1_hac.cfg \
--compress_fastq

cat ./basecalled/pass/*.fastq.gz > "./basecalled/${indir_basename}.fastq.gz"

# After the job is done we copy our output back to our directory
outdir="/home/Magda_lab/maitenat/scratch60/ont/${indir_basename}"
mkdir -p ${outdir} && cp -r "${SCRATCH_DIRECTORY}/basecalled/${indir_basename}.fastq.gz" ${outdir}
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
