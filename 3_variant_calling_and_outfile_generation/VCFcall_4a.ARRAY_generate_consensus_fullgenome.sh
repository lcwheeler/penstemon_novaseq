#!/bin/sh

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p defq
#SBATCH --job-name=variants_bcftools

#wessinger-48core

###NOTE: this is set up as a batch array script. SLURM submit array code begins on line 47.


cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

########################
#source appropriate environments to enable use of conda installs through batch submission
source /home/lw74/.bashrc
source /home/lw74/.bash_profile

#activate conda environment with packages installed
#needs samtools (v1.15.1 works) and bcftools (v1.15.1 works)
#conda activate mapping_etc


#hard path to the reference genome and the samtools faidx index of the genome
#mapped filtered reads should not have a "/" at the end. This should be the same directory as output from script: mapping_pipeline_array_s1a-f.sh
refgenome="/work/lw74/refGenome/Pbar.2022.LG.fa"
faidxfile="/work/lw74/refGenome/Pbar.2022.LG.fa.fai"
mapped_filtered_reads="/work/lw74/habro/mapped_filtered_bams"


#Location of filtered vcf, from script #3 in this pipeline.
invcf="/work/lw74/habro/VCFs_WithOutgroups/filtered_consensus_ready_no-indels.vcf.gz"


#path to the output directory for consensus sequence alignments. Make prior to running.
mkdir -p "/work/lw74/habro/consensus_alignments_WithOutgroups/individual_fullgenome_fastas"
outdir="/work/lw74/habro/consensus_alignments_WithOutgroups/individual_fullgenome_fastas"


#number of cores used
numthreads=8


#set up list of samples for which to generate consensus sequences
#this effectively lists sample names as listed in the .vcf, if following same protocol
#if sample names were changed at some point in .vcf, this will throw an error
#infile is the full path to each input bam file
#insample is just the sample name (no path)
#for batch array, sbatch --array [0-n], where n is number of .bam files in $mapped_filtered_reads -1
filelist=($mapped_filtered_reads/*.bam)
infile="${filelist[$SLURM_ARRAY_TASK_ID]}"

insample=${infile%_mapped_filtered.bam}
insample=${insample##*/}

#get list of scaffold names from samtools index
scaffoldnames=$(awk '{print $1}' $faidxfile)

########################
#1. index the reference for annotations of interest (just scaffold names if whole genome)
#2. pipe this into bcftools consensus to generate sequence alignment for samples


#generate consensus sequence
#--haplotype I outputs IUPAC ambiguity codes based on sample genotypes
#--absent changes all sites not included in the vcf to "N" (rather than to consensus)
#--missing changes all ./. sites in vcf to "N"
#--mark-del changes deletions to "-"
#--regions specifies which scaffolds you want to consider (if you don't specify, it only does one scaffold)
#--sample is specifying which sample to generate sequence for -- this is why array is used

samtools faidx $refgenome ${scaffoldnames[@]} | 
 bcftools consensus --haplotype I \
 --absent "N" \
 --missing "N" \
 --mark-del "-" \
 --sample $infile \
 $invcf > $outdir/consensus_fullgenome_$insample.fa


