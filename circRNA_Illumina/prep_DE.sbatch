#!/bin/bash
#SBATCH --job-name=ciriDE2
#SBATCH --output=ciriDE2_%A.out
#SBATCH --error=ciriDE2_%A.err
#SBATCH --ntasks=1
#SBATCH --partition=q_cn
#SBATCH --mem-per-cpu=40G

# Description:  Prepare stringtie output to afterwards do a differential expression analysis of circRNAs using biological replicates
# It is the second step after preparing CIRIquant output

module load stringtie/2.1.2

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}

sample_list=$1 # this is the _out list file
sample_list_basename=`basename $sample_list`
exp_name="${sample_list_basename%%.*}"

prepDE.py -i $sample_list

# After the job is done we copy our output back to our directory
cp gene_count_matrix.csv ${exp_name}_gene_count_matrix.csv
outdir="/GPFS/Magda_lab_temp/maitenat/illumina_circ/ciriquant"
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/${exp_name}_gene_count_matrix.csv ${outdir}

# Clean up after ourselves
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}