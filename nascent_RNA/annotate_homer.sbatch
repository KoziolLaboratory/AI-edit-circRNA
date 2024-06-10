#!/bin/bash
#SBATCH --job-name=homerann
#SBATCH --output=annotate_homer_%A.out
#SBATCH --error=annotate_homer_%A.err
#SBATCH --ntasks=1
#SBATCH --partition=q_cn
#SBATCH --mem-per-cpu=20G

# Description:  Annotate edits into their genomic locations. 

module load homer/homer

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
filename=$1

# Copy necessary files to scratch directory
cp $filename ${SCRATCH_DIRECTORY}
file_basename=`basename $filename`
samplename="${file_basename%%.*}"
# Genome
cp /home/Magda_lab/maitenat/DATA/indeces/mm* ${SCRATCH_DIRECTORY}
cp /home/Magda_lab/maitenat/DATA/mm39.chromFa.tar.gz ${SCRATCH_DIRECTORY}
tar -xzvf mm39.chromFa.tar.gz
cat chr*.fa > mm39.fa
rm chr*
cp /GPFS/Magda_lab_permanent/maitenat/mm39.ncbiRefSeq.gtf ${SCRATCH_DIRECTORY}

# Convert to bed file
tail -n +2 ${file_basename} > ${samplename}.temp0
cat ${samplename}.temp0 | awk -v OFS='\t' '{ print $1, ($2-1), $2, $1"_"$2, $3, $4 }' > ${samplename}.temp1
awk -v OFS='\t' '{if ($6=="0") {$6="-"} else {$6="+"}}1' ${samplename}.temp1 > ${samplename}.temp
mv ${samplename}.temp ${samplename}.bed

# Run homer
parseGTF.pl mm39.ncbiRefSeq.gtf ann -annTSSstartOffset -1 -annTSSendOffset 1 -annTTSstartOffset -1 -annTTSendOffset 1 > annotations.txt
cat annotations.txt | awk '$1 == "exon"' > annotations_sorted.txt
cat annotations.txt | awk '$1 == "intron"' >> annotations_sorted.txt
cat annotations.txt | awk '$1 == "promoter-TSS"' >> annotations_sorted.txt
cat annotations.txt | awk '$1 == "TTS"' >> annotations_sorted.txt

assignGenomeAnnotation annotations_sorted.txt annotations_sorted.txt -prioritize annotations.final.txt > stats.txt
annotatePeaks.pl ${samplename}.bed mm39.fa -ann annotations.final.txt > ${samplename}_annot.bed


# After the job is done we copy our output back to our directory
outdir="/home/Magda_lab/maitenat/scratch60/nascent_RNA/reditools/"
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/${samplename}_annot.bed ${outdir}
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
