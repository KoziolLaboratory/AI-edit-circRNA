#!/bin/bash
#SBATCH --job-name=ciri2
#SBATCH --output=ciri2_%A.out
#SBATCH --error=ciri2_%A.err
#SBATCH --ntasks=1
#SBATCH --partition=q_fat
#SBATCH --mem-per-cpu=180G

# Description: Run CIRI2 on a sam file to detect circRNAs

module load ciri/2.0.6

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
samfile=$1

# Copy necessary files to scratch directory
cp $samfile ${SCRATCH_DIRECTORY}
samfile_basename=`basename $samfile`
echo ${samfile_basename}
samplename="${samfile_basename%%.*}"

cp /GPFS/Magda_lab_permanent/maitenat/DATA/mm39.ncbiRefSeq.gtf ${SCRATCH_DIRECTORY}
cp /home/Magda_lab/maitenat/DATA/indeces/mm* ${SCRATCH_DIRECTORY}
cp /home/Magda_lab/maitenat/DATA/mm39.chromFa.tar.gz ${SCRATCH_DIRECTORY}
tar -xzvf mm39.chromFa.tar.gz
echo "Uncompressed!"
cat chr*.fa > mm39.fa
rm chr*
rm mm39.chromFa.tar.gz

# Run CIRI2
perl $CIRI_HOME/CIRI2.pl -I ${samfile_basename} -O ${samplename}_ciri2.txt -F mm39.fa -A mm39.ncbiRefSeq.gtf

# After the job is done we copy our output back to our directory
outdir="/GPFS/Magda_lab_temp/maitenat/illumina_circ/ciri2"
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/${samplename}_ciri2.txt ${outdir}

cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}