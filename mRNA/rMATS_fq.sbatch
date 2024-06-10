#!/bin/bash
#SBATCH --job-name=rMATS
#SBATCH --output=rmats_%A.out
#SBATCH --error=rmats_%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=q_cn
#SBATCH --mem-per-cpu=20G

# Description:  Run rMATS to detect differential splicing events

module load rMATS/4.1.2

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
b1=/home/Magda_lab/maitenat/scratch60/mRNA_illumina/b1_fq.txt
b2=/home/Magda_lab/maitenat/scratch60/mRNA_illumina/b2_fq.txt
gtf_file=/GPFS/Magda_lab_permanent/maitenat/mm39.ncbiRefSeq.gtf
star_idx=/GPFS/Magda_lab_permanent/maitenat/indeces/STAR

# Run rMATS
mkdir -p tmp_output
mkdir -p output
python $rMATS/rmats.py --s1 ${b1} \
--s2 ${b2} \
--bi ${star_idx} \
--gtf ${gtf_file} \
--readLength 150 \
--nthread 4 \
--allow-clipping \
--od ./output \
--tmp ./tmp_output

# After the job is done we copy our output back to our directory
outdir="/home/Magda_lab/maitenat/scratch60/mRNA_illumina/splicing"
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/output/* ${outdir}

cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
