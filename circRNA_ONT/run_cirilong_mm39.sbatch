#!/bin/bash
#SBATCH --job-name=cirilong
#SBATCH --output=cirilong_%A.out
#SBATCH --error=cirilong_%A.err
#SBATCH --ntasks=1
#SBATCH --partition=q_fat
#SBATCH --mem-per-cpu=80G

# Description:  This script runs CIRI-long to identify circular RNAs from ONT sequencing reads

module load ciri-long/1.0.3
module load bwa/0.7.17
module load samtools/1.12

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
infile=$1 # fastq file

# Copy necessary files to scratch directory
cp $infile ${SCRATCH_DIRECTORY}
infile_basename=`basename $infile`
echo ${infile_basename}
samplename="${infile_basename%%.*}"

cp /GPFS/Magda_lab_permanent/maitenat/DATA/mm39.ncbiRefSeq.gtf ${SCRATCH_DIRECTORY}
cp /home/Magda_lab/maitenat/DATA/indeces/* ${SCRATCH_DIRECTORY}
cp /home/Magda_lab/maitenat/DATA/mm39.chromFa.tar.gz ${SCRATCH_DIRECTORY}
tar -xzvf mm39.chromFa.tar.gz
echo "Uncompressed!"
cat chr*fa > mm39.fa
rm chr*
rm mm39.chromFa.tar.gz

# Filter fastq file by readlength
zcat ${infile_basename} | awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 500) {print header, seq, qheader, qseq}}' - > temp.fastq
echo "Number of filtered reads (remaining):"
cat temp.fastq | echo $((`wc -l`/4))
mv temp.fastq ${samplename}.fastq

# Run CIRI-long to identify circular reads from sequencing reads
CIRI-long call -i ${samplename}.fastq \
               -o ./${samplename}_call \
               -r mm39.fa \
               -p ${samplename} \
               -t 4 \
               -a mm39.ncbiRefSeq.gtf

touch $samplename.lst
echo "${samplename} ./${samplename}_call/${samplename}.cand_circ.fa" >> ${samplename}.lst

# Change the names to those contigs with "-" in their name
sed -i 's/SVEVBRD-HPRT-B-M2/SVEVBRD_HPRT_B_M2/' ./${samplename}_call/${samplename}.cand_circ.fa
sed -i 's/SVEVBRD-HPRT-B-M2/SVEVBRD_HPRT_B_M2/' ./${samplename}_call/${samplename}.low_confidence.fa
# Need to do it also for the reference genome and indeces
sed -i 's/SVEVBRD-HPRT-B-M2/SVEVBRD_HPRT_B_M2/' mm39.fa
sed -i 's/SVEVBRD-HPRT-B-M2/SVEVBRD_HPRT_B_M2/' mm39.fa.amb
sed -i 's/SVEVBRD-HPRT-B-M2/SVEVBRD_HPRT_B_M2/' mm39.fa.ann
sed -i 's/SVEVBRD-HPRT-B-M2/SVEVBRD_HPRT_B_M2/' mm39.fa.bwt
sed -i 's/SVEVBRD-HPRT-B-M2/SVEVBRD_HPRT_B_M2/' mm39.fa.fai
sed -i 's/SVEVBRD-HPRT-B-M2/SVEVBRD_HPRT_B_M2/' mm39.fa.pac
sed -i 's/SVEVBRD-HPRT-B-M2/SVEVBRD_HPRT_B_M2/' mm39.fa.sa

CIRI-long collapse -i ./${samplename}.lst \
                    -o ./${samplename}_collapse \
                    -p ${samplename} \
                    -r mm39.fa \
                    -a mm39.ncbiRefSeq.gtf

# After the job is done we copy our output back to our directory
outdir=/home/Magda_lab/maitenat/scratch60/ont/cirilong_mm39
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/${samplename}_call ${outdir}
cp -r ${SCRATCH_DIRECTORY}/${samplename}_collapse ${outdir}
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
