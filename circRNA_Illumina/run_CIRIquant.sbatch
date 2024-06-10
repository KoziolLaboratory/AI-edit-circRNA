#!/bin/bash
#SBATCH --job-name=ciriquant
#SBATCH --output=ciri-quant_%A.out
#SBATCH --error=ciri-quant_%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=q_fat
#SBATCH --mem-per-cpu=80G

# Description:  Run ciri-quant to quantify circRNAs in RNA-Seq data
# This is done per sample independently

module load ciriquant/1.1.2
module load bwa/0.7.17
module load hisat2/2.1.0
module load stringtie/2.1.2
module load samtools/1.12

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}

fq_1=$1
fq_2=$2
ciri2_results=$3
config=/GPFS/Magda_lab_temp/maitenat/illumina_circ/config.yml

# Copy necessary files to scratch directory
cp $fq_1 ${SCRATCH_DIRECTORY}
fq_1_basename=`basename $fq_1`
cp $fq_2 ${SCRATCH_DIRECTORY}
fq_2_basename=`basename $fq_2`
cp $ciri2_results ${SCRATCH_DIRECTORY}
ciri2_results_basename=`basename $ciri2_results`
samplename="${fq_1_basename%%.*}"

# Run CIRIquant
CIRIquant -t 4 \
          -1 ${fq_1_basename} \
          -2 ${fq_2_basename} \
          --config $config \
          -o ./${samplename} \
          -p $samplename \
          --circ $ciri2_results \
          --tool CIRI2


# After the job is done we copy our output back to our directory
outdir="/GPFS/Magda_lab_temp/maitenat/illumina_circ/ciriquant"
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/${samplename}/*.gtf ${outdir}
cp -r ${SCRATCH_DIRECTORY}/${samplename}/gene/*list ${outdir}
cp -r ${SCRATCH_DIRECTORY}/${samplename}/gene/*out.gtf ${outdir}

# Clean up after ourselves
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
