#!/bin/sh

#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p defq
#SBATCH --job-name=filter_invariants

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
#need bcftools (v1.15.1 works)
#conda activate mapping_etc

#load vcftools (v0.1.17 on cluster works)
module load vcftools


#Location of unfiltered vcf, from script #2 in this pipeline.
invcf="/work/lw74/habro/VCFs_WithOutgroups/unfiltered_vcf.gz"


#path to the output directory for vcfs and gvcfs. Should be made already.
outdir="/work/lw74/habro/VCFs_WithOutgroups"


#number of cores used converting from gvcf to vcf
numthreads=16


########################

#####
#VCF FILTERING#
#####

#We must filter variant and invariant sites separately -- applying the same filters to all sites would remove invariant sites because of the way they are coded in the file.


#only invariant sites.
#max-maf 0 means no sites with minor allele frequency > 0.
#we will not implement any filters here (but you could feasibly do min depth filters)
#pipe into view because the compression is much faster than bgzip
 vcftools --gzvcf $invcf \
 --max-maf 0 \
 --recode --recode-INFO-all --stdout | 
 bcftools view - \
 --threads $numthreads \
 -Oz -o $outdir/tmp-invariants.filtered.vcf.gz


