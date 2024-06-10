#!/bin/bash
#SBATCH --job-name=ciriDE3
#SBATCH --output=ciriDE3_%A.out
#SBATCH --error=ciriDE3_%A.err
#SBATCH --ntasks=1
#SBATCH --partition=q_cn
#SBATCH --mem-per-cpu=40G

# Description:  Perform differential expression analysis of circRNAs using biological replicates
# It is the third step after preparing CIRIquant output and stringtie output

module load R/4.1.0
module load ciriquant/1.1.2

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}

library=$1
bsj=$2
gene=$3
library_basename=`basename $library`
exp_name=${library_basename//_library_info.csv/}

CIRI_DE_replicate --lib  $library \
            --bsj  $bsj \
            --gene $gene \
            --out  ${exp_name}_DE.tsv

# After the job is done we copy our output back to our directory
outdir="/GPFS/Magda_lab_temp/maitenat/illumina_circ/ciriquant"
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/${exp_name}_DE.tsv ${outdir}

# Clean up after ourselves
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
