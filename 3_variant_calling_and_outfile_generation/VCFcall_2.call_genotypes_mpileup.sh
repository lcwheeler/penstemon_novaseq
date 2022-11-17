#!/bin/sh

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p defq
#SBATCH --job-name=variants_bcftools

#wessinger-48core
cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

########################
#source appropriate environments to enable use of conda installs through batch submission
source /home/lw74/.bashrc
source /home/lw74/.bash_profile

#activate conda environment with packages installed
#needs bcftools (v1.15.1 works)
#conda activate mapping_etc


#hard path to the reference genome and mapped, filtered reads
#mapped filtered reads should be same directory as output from script: mapping_pipeline_array_s1a-f.sh
refgenome="/work/lw74/refGenome/Pbar.2022.LG.fa"
mapped_filtered_reads="/work/lw74/habro/mapped_filtered_bams"


#path to the output directory for vcf files (called variants). Make prior to running
mkdir -p "/work/lw74/habro/VCFs_WithOutgroups"
outdir="/work/lw74/habro/VCFs_WithOutgroups"


#number of cores used for mpileup and call
#change this to match number of cores allocated
numthreads=8
########################


#make list of mapped, filtered bams
bamlist=($mapped_filtered_reads/*.bam)


#Produce GT likelihoods, call variants, and normalize indels -> unfiltered .vcf
#mpileup produces genotype likelihoods from bam files. Filtering for BQ and MQ > 20
#--gvcf call contiguous reference haplotype blocks. require at least 2 reads for inclusion
#call calls SNPs from the genotype likelihoods
bcftools mpileup --threads $numthreads -Ou \
 -a FORMAT/AD,FORMAT/DP \
 --skip-indels \
 --min-BQ 20 \
 --min-MQ 20 \
 -f $refgenome ${bamlist[@]} | 
 bcftools call --threads $numthreads -m --gvcf 2 \
 -Oz -o $outdir/unfiltered_vcf.gz


#index the vcf
tabix $outdir/unfiltered_vcf.gz

