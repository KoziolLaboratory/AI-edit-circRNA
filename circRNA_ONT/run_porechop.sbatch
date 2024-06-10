#!/bin/bash
#SBATCH --job-name=porechop
#SBATCH --output=porechop_%A.out
#SBATCH --error=porechop_%A.err
#SBATCH -n 1
#SBATCH --partition=q_cn
#SBATCH --mem-per-cpu=30G

# Description: Remove barcodes and adapters from Nanopore reads

module load porechop/0.2.4
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
samplename="${infile_basename%%.*}"

# Remove barcodes and adapters
porechop -i $infile_basename -o $samplename"_clean.fastq.gz"
echo "porechop finished"

# Run FastQC
mkdir fastqc
fastqc -t 10 -o fastqc $samplename"_clean.fastq.gz"

# After the job is done we copy our output back to our directory
outdir=/home/Magda_lab/maitenat/scratch60/ont
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/*_clean.fastq.gz ${outdir}
cp -r ${SCRATCH_DIRECTORY}/fastqc/* ${outdir}
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
exit 0
