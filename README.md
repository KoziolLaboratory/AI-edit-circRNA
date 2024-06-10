# A-to-I editing regulates splicing and circRNAs in the brain

This repository contains the code associated with the article titled "A-to-I editing regulates splicing and circRNAs in the brain". The code provided here was used for the analysis and generation of results presented in the article.

## Introduction

The code in this repository is organized into four main folders:

1. `circRNA_Illumina`: This folder contains the code used for the identification of circRNAs in Illumina RNA sequencing data.

2. `circRNA_ONT`: This folder contains the code used for the identification of circRNAs in Oxford Nanopore Technologies (ONT) sequencing data.

3. `mRNA`: This folder contains scripts related to Illumina mRNA sequencing data.

4. `nascent_RNA`: This folder contains scripts related to the analysis of nascent RNA sequencing data.


## Code Overview and usage

The scripts used in this work are written in Shell (`.sbatch`), R (`.R`), and possibly others depending on the specific requirements of the analysis. The `.sbatch` files are scripts written for SLURM job scheduler for Linux. These scripts are used to submit jobs to the SLURM scheduler, which then manages the job's execution on the cluster. 

### circRNA_Illumina

The `circRNA_Illumina` folder contains scripts for the identification of circRNAs in Illumina RNA sequencing data. The starting data for these scripts are clean fastq files. The scripts in this section should be run in the order listed.

1. `run_bwa.sbatch`
2. `hisat_index.sbatch`
3. `run_CIRI2.sbatch`
4. `run_CIRIquant.sbatch`
5. `prep_CIRIquant.sbatch`
6. `prep_DE.sbatch`
7. `run_CIRI_DE.sbatch`


### circRNA_ONT

The `circRNA_ONT` folder contains the data for the identification of circRNAs in ONT sequencing data. The starting data for these scripts are fast5 files generated by ONT sequencing chips. The scripts in this section should be run in the order listed.

1. `run_basecalling.sbatch`
2. `run_fastqc.sbatch`
3. `run_porechop.sbatch`
4. `run_cirilong_mm39.Rmd`


### mRNA

The `mRNA` folder contains a few scripts for splicing analyses on mRNA Illumina sequencing data. For the analysis of this data, the following scripts are called, in this order:

1. `run_trimmomatic.sbatch`
2. `rMATS_fq.sbatch`

### nascent_RNA

The `nascent_RNA` folder contains scripts to identify A-to-I edits in nascent RNA sequencing data. These scripts need to be called in this order:

1. `run_fastqc.sbatch`
2. `run_trimmomatic.sbatch`
3. `STAR_align.sbatch`
4. `run_reditools.sbatch`
5. `filter_reditools.sbatch`
6. `annotate_homer.sbatch`


## Requirements

This project requires the following tools and libraries:

### Software

- [R](https://www.r-project.org/)
- [guppy](https://community.nanoporetech.com/protocols/Guppy-protocol/v/GPB_2003_v1_revn_14Dec2018) for basecalling
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for quality control checks on raw sequence data
- [porechop](https://github.com/rrwick/Porechop) for adapter trimming in ONT data
- [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) for read trimming
- [HISAT2](http://daehwankimlab.github.io/hisat2/) for read alignment
- [BWA](http://bio-bwa.sourceforge.net/) for read alignment
- [stringtie](https://ccb.jhu.edu/software/stringtie/) for transcript assembly and quantification
- [samtools](http://www.htslib.org/) for manipulating aligned reads
- [rMATS](http://rnaseq-mats.sourceforge.net/) for differential alternative splicing analysis
- [CIRIquant](https://ciri-cookbook.readthedocs.io/en/latest/CIRIquant_0_home.html) for circRNA quantification
- [CIRI2](https://ciri-cookbook.readthedocs.io/en/latest/CIRI2.html) for circRNA detection
- [CIRI-long](https://ciri-cookbook.readthedocs.io/en/latest/CIRI-long_0_home.html) for long-read circRNA detection
- [homer](http://homer.ucsd.edu/homer/) for variant annotation
- [reditools](https://sourceforge.net/projects/reditools/) for RNA editing detection
- [bcftools](http://www.htslib.org/doc/bcftools.html) for manipulating variant calls


### Data
- Reference genome (GRCm39/mm39)
- Adapter sequences for the sequencing platform used
- Gene annotations (NCBI RefSeq: mm39.ncbiRefSeq.gtf)
- Variant annotations (mgp_REL2021_snps.bed, GCA_000001635.9_current_ids.bed, mm10-blacklist.v2.Liftover.mm39.bed.txt)


### Job Scheduler
- [SLURM](https://slurm.schedmd.com/overview.html) for job scheduling on a high-performance computing cluster. Many of the scripts in this project are written as SLURM batch scripts and require SLURM for execution. However, these SLURM scripts are essentially shell scripts. If you are not using SLURM, you can adapt them to run as regular bash scripts with some minimal knowledge of bash scripting.

Please ensure all dependencies are installed and the raw data is available in the appropriate directories before running the scripts.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
