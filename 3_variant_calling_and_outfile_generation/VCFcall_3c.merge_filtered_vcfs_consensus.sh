#!/bin/sh

#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p defq
#SBATCH --job-name=merge_vcf_consensus
#SBATCH --output merge_vcf_consensus%j.out
#SBATCH --error merge_vcf_consensus%j.err


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


#index both vcfs
tabix $outdir/tmp-invariants.filtered.vcf.gz
tabix $outdir/tmp-variants.filtered.vcf.gz


#concatenate the two vcfs and remove tmp files
bcftools concat $outdir/tmp-invariants.filtered.vcf.gz \
 $outdir/tmp-variants.filtered.vcf.gz \
 --threads $numthreads --allow-overlaps -Oz \
 -o $outdir/filtered_consensus_ready_no-indels.vcf.gz

#you could uncomment this to remove the tmp files, but I prefer to remove by hand, in case something goes wrong with the script.
#rm $outdir/tmp-invariants.filtered.bcf.gz*
#rm $outdir/tmp-variants.filtered.bcf.gz*


#index final merged file
tabix $outdir/filtered_consensus_ready_no-indels.vcf.gz

