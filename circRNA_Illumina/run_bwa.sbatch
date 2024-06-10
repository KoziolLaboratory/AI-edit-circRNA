#!/bin/bash
#SBATCH --job-name=align
#SBATCH --output=align_%A.out
#SBATCH --error=align_%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=q_fat
#SBATCH --mem-per-cpu=40G

# Description:  Aligns Illumina reads to the mouse genome using BWA

# We are keeping SAM files as required by CIRI2

module load bwa/0.7.17
module load samtools/1.12

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
R1=$1
R2=$2

# Copy necessary files to scratch directory
cp $R1 ${SCRATCH_DIRECTORY}
cp $R2 ${SCRATCH_DIRECTORY}
R1_basename=`basename $R1`
R2_basename=`basename $R2`
sample="${R1_basename%%.*}"

cp /home/Magda_lab/maitenat/DATA/indeces/mm* ${SCRATCH_DIRECTORY}
cp /home/Magda_lab/maitenat/DATA/mm39.chromFa.tar.gz ${SCRATCH_DIRECTORY}
tar -xzvf mm39.chromFa.tar.gz
echo "Uncompressed!"
cat chr*.fa > mm39.fa
rm chr*
rm mm39.chromFa.tar.gz

# Run alignment
echo "Running alignment..."
bwa mem -t 4 mm39.fa ${R1_basename} ${R2_basename} > ${sample}.sam

# After the job is done we copy our output back to our directory
outdir="/GPFS/Magda_lab_temp/maitenat/illumina_circ/"
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/${sample}.sam ${outdir}

cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
