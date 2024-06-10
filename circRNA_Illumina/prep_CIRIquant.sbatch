#!/bin/bash
#SBATCH --job-name=ciriDE1
#SBATCH --output=ciriDE1_%A.out
#SBATCH --error=ciriDE1_%A.err
#SBATCH --ntasks=1
#SBATCH --partition=q_cn
#SBATCH --mem-per-cpu=40G

# Description:  Prepare CIRIquant output to afterwards do a differential expression analysis of circRNAs using biological replicates

module load ciriquant/1.1.2
module load samtools/1.12

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}

sample_list=$1
sample_list_basename=`basename $sample_list`
exp_name="${sample_list_basename%%.*}"

prep_CIRIquant -i $sample_list \
                --lib ${exp_name}_library_info.csv \
                --circ ${exp_name}_circRNA_info.csv \
                --bsj ${exp_name}_bsj.csv \
                --ratio ${exp_name}_ratio.csv

# After the job is done we copy our output back to our directory
outdir="/GPFS/Magda_lab_temp/maitenat/illumina_circ/ciriquant"
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/*.csv ${outdir}

# Clean up after ourselves
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
