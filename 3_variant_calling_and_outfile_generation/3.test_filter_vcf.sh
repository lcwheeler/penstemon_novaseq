#!/bin/sh

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p wessinger-48core
#SBATCH --job-name=filter_vcf_consensus-ready


cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

########################
#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile


#activate conda environment with packages installed
#needs samtools (v1.15.1 works) and bcftools (v1.15.1 works)
conda activate mapping_etc

#load vcftools (v0.1.17 on cluster works)
module load vcftools


#Location of unfiltered vcf, from script #2 in this pipeline.
invcf="/work/bs66/dasanthera_novaseq/GVCF_VCFs/unfiltered_vcf.gz"


#path to the output directory for vcfs and gvcfs. Should be made already.
outdir="/work/bs66/dasanthera_novaseq/GVCF_VCFs"


#number of cores used converting from gvcf to vcf
numthreads=8


########################

#####
#VCF FILTERING#
#####

#We must filter variant and invariant sites separately -- applying the same filters to all sites would remove invariant sites because of the way they are coded in the file.


#only invariant sites.
#max-maf 0 means no sites with minor allele frequency > 0.
#we don't need filters on these sites -- just pull them out so they aren't removed when filtering the variant sites.
 vcftools --gzvcf $invcf \
 --max-maf 0 \
 --recode --recode-INFO-all --stdout | 
 bgzip -c > $outdir/tmp-invariants_filtered.vcf.gz


#only variant sites. filters on QUALITY, DEPTH, and MISSIGNESS of reads.
#note that exact parameter values should be thought out on analysis-by-analysis basis. This script is for filtering with the eventual goal of generating multiple sequence alignments (phylogenomics!). So, some of the population-genomic filters common at the vcf filtering stage are not helpful here (e.g., filtering on minor allele frequency, HWE proportions, etc.).

#--min-alleles 2 means there must be at least two alleles (variant)
#--minGQ recode genotypes below GQ [int] as missing
# meanDP filters remove sites with mean allelic depths outside specified limits [int]
#--max-missing removes sites on basis of missing data (1 = no missing data allowed)
 vcftools --gzvcf $invcf \
 --min-alleles 2 \
 --minGQ 15 \
 --min-meanDP 3 \
 --max-meanDP 30 \
 --max-missing 0.80 \
 --recode --recode-INFO-all --stdout | 
 bgzip -c > $outdir/tmp-variants_filtered.vcf.gz


#index both vcfs
tabix $outdir/tmp-invariants_filtered.vcf.gz
tabix $outdir/tmp-variants_filtered.vcf.gz


#concatenate the two vcfs and remove tmp files
bcftools concat $outdir/tmp-invariants_filtered.vcf.gz \
 $outdir/tmp-variants_filtered.vcf.gz \
 --threads 4 -a -Oz \
 -o $outdir/filtered_consensus-ready.vcf.gz

rm $outdir/tmp-invariants_filtered.vcf.gz*
rm $outdir/tmp-variants_filtered.vcf.gz*


#index final merged file
tabix $outdir/filtered_consensus-ready.vcf.gz



