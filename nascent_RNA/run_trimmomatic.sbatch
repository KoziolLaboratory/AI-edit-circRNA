#!/bin/bash
#SBATCH --job-name=trimmo
#SBATCH --output=trimmomatic_%A.out
#SBATCH --error=trimmomatic_%A.err
#SBATCH --ntasks=1
#SBATCH --partition=q_fat
#SBATCH --mem-per-cpu=120G

# Description: Trimming of adapters and low quality reads using trimmomatic

module load trimmomatic/0.39

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

# Run trimmomatic
echo "Trimming adapters on sample ${sample}"
trimmomatic PE -threads 6 $R1 $R2 $sample"_trimmed_R1.fastq" $sample"_unpaired_R1.fastq" $sample"_trimmed_R2.fastq" $sample"_unpaired_R2.fastq" ILLUMINACLIP:/usr/nzx-cluster/apps/trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10:10:TRUE SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:36
gzip $sample"_trimmed_R1.fastq"
gzip $sample"_trimmed_R2.fastq"

# After the job is done we copy our output back to our directory
outdir="/home/Magda_lab/maitenat/scratch60/nascent_RNA/clean"
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/*trimmed*.fastq.gz ${outdir}

cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
