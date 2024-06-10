#!/bin/bash
#SBATCH --job-name=STAR
#SBATCH --output=STAR-align_%A.out
#SBATCH --error=STAR-align_%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=q_fat
#SBATCH --mem=140G

# Description: Aligns reads to the genome using STAR

module load rMATS/4.1.2
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
sample=`echo "${R1_basename//_trimmed_R1.fastq.gz/}"`

mkdir -p STAR
cp /home/Magda_lab/maitenat/DATA/indeces/STAR/* ${SCRATCH_DIRECTORY}/STAR


# Run alignment
echo "Running alignment..."
/usr/nzx-cluster/apps/rMATS/4.1.2/rmats/bin/STAR --runThreadN 4 \
--genomeDir ./STAR \
--genomeLoad NoSharedMemory \
--outFileNamePrefix ${sample}_ \
--outSAMtype BAM SortedByCoordinate \
--readFilesCommand zcat \
--readFilesIn ${R1_basename} ${R2_basename}

# Index aligned file
samtools index ${sample}_Aligned.sortedByCoord.out.bam

# After the job is done we copy our output back to our directory
outdir="/home/Magda_lab/maitenat/scratch60/nascent_RNA/aligned"
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/${sample}_Aligned.sortedByCoord.out.ba* ${outdir}

cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
