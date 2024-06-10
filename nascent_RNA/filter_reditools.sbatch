#!/bin/bash
#SBATCH --job-name=filtREDI
#SBATCH --output=filterREDI_%A.out
#SBATCH --error=filterREDI_%A.err
#SBATCH --ntasks=1
#SBATCH --partition=q_fat
#SBATCH --mem-per-cpu=80G

# Description: Filter REDItools output to keep only A to I editing sites with coverage > 10 and not in blacklist

module load R/3.6.2
module load bcftools/1.15.1

# Parameters of the current job
pwd; hostname; date

# Define and create a unique scratch directory for this job
SCRATCH_DIRECTORY=/var/tmp/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
varfile=$1

R_script=/GPFS/Magda_lab_permanent/maitenat/Illumina_RNASeq/code/annotate_SNPs.R
dbsnp=/GPFS/Magda_lab_permanent/maitenat/mgp_REL2021_snps.bed
blacklist=/GPFS/Magda_lab_permanent/maitenat/mm10-blacklist.v2.Liftover.mm39.bed.txt
eva_snp=/GPFS/zhangli_lab_permanent/zhouxiangyu/novo/CIBR-MAGDA/20221213/GCA_000001635.9_current_ids.vcf.gz
awk_script=/GPFS/Magda_lab_permanent/maitenat/Illumina_RNASeq/code/filter_pos_range.awk
cp ${varfile} ${SCRATCH_DIRECTORY}

var_basename=`basename $varfile`
sample="${var_basename%%.*}"
sample_cut="${sample%_*}"

# Filter A to I editing sites
head -n 1 ${var_basename} > header
cat header > ${sample_cut}_AI.txt
awk '$11=="AG"' ${var_basename} >> ${sample_cut}_AI.txt

# Filter coverage
awk '$5>=10' ${sample_cut}_AI.txt > ${sample_cut}_AI_cov.txt

# Filter SNPs
## dbSNP
awk '{print $1":"$3}' ${dbsnp} > dbsnp_id
## EVA
bcftools view -H ${eva_snp} > eva_clean.vcf
cat eva_clean.vcf | awk 'length($4) == 1 && length($5) == 1' > eva_clean0.vcf
mv eva_clean0.vcf eva_clean.vcf
awk '{print "chr"$1":"$2}' eva_clean.vcf >> allsnp_id

R --no-echo --no-restore --file=${R_script} --args ${sample_cut}_AI_cov.txt allsnp_id

# Filter blacklist
awk -f ${awk_script} ${blacklist} ${sample_cut}_AI_cov_noSNP.txt > muts_in_blacklist
sort muts_in_blacklist > muts_in_blacklist_sorted
sort ${sample_cut}_AI_cov_noSNP.txt > ${sample_cut}_AI_cov_noSNP_sorted.txt
# Get the difference between the two files
comm -23 ${sample_cut}_AI_cov_noSNP_sorted.txt muts_in_blacklist_sorted > ${sample_cut}_AI_cov_noSNP_noblack.txt

# After the job is done we copy our output back to our directory
outdir="/home/Magda_lab/maitenat/scratch60/nascent_RNA/reditools/"
mkdir -p ${outdir} && cp -r ${SCRATCH_DIRECTORY}/${sample_cut}_AI_cov_noSNP_noblack.txt ${outdir}
cp -r ${SCRATCH_DIRECTORY}/${sample_cut}_AI_cov_noSNP.txt ${outdir}

cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}
