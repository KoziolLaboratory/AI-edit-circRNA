#!/bin/bash
#SBATCH --job-name=reditool
#SBATCH --output=reditools_%A.out
#SBATCH --error=reditools_%A.err
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --partition=q_fat
#SBATCH --mem-per-cpu=80G

# Description: Run reditools to detect RNA editing sites

module load reditools/2.0

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
bamfile=$1

# Copy necessary files to scratch directory
cp $bamfile ${SCRATCH_DIRECTORY}
cp ${bamfile}".bai" ${SCRATCH_DIRECTORY}
bam_basename=`basename $bamfile`
sample="${bam_basename%%.*}"

cp /home/Magda_lab/maitenat/DATA/indeces/* ${SCRATCH_DIRECTORY}
cp /home/Magda_lab/maitenat/DATA/mm39.chromFa.tar.gz ${SCRATCH_DIRECTORY}
tar -xzvf mm39.chromFa.tar.gz
echo "Uncompressed!"
cat chr* > mm39.fa
rm chr*
rm mm39.chromFa.tar.gz

# Run reditools
python /usr/nzx-cluster/apps/reditools2/2.0/reditools2.0/src/cineca/reditools.py -f ${bam_basename} \
-r mm39.fa \
-c \
-os 5 \
-s 2 \
-q 10 \
-bq 20 \
-mrl 50 \
-m homop.bed \
-o ${sample}_stats.txt \
--strict \
-men 3 \
--max-editing-nucleotides 2 \
-C

# After the job is done we copy our output back to our directory
outdir="/home/Magda_lab/maitenat/scratch60/nascent_RNA/reditools"
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/${sample}_stats.txt ${outdir}
cp -r homop.bed ${outdir}

cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}

