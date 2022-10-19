#!/bin/sh

#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p wessinger-48core
#SBATCH --job-name=merge_vcf_consensus


cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

########################
#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile


#activate conda environment with packages installed
#need bcftools (v1.15.1 works)
conda activate mapping_etc

#load vcftools (v0.1.17 on cluster works)
module load vcftools


#Location of unfiltered vcf, from script #2 in this pipeline.
invcf="/work/bs66/dasanthera_novaseq/GVCF_VCFs/unfiltered_vcf.gz"


#path to the output directory for vcfs and gvcfs. Should be made already.
outdir="/work/bs66/dasanthera_novaseq/GVCF_VCFs"


#number of cores used converting from gvcf to vcf
numthreads=16


########################

#####
#VCF FILTERING#
#####


#index both vcfs
tabix $outdir/tmp-invariants.filtered.bcf.gz
tabix $outdir/tmp-variants.filtered.bcf.gz


#concatenate the two vcfs and remove tmp files
bcftools concat $outdir/tmp-invariants.filtered.bcf.gz \
 $outdir/tmp-variants.filtered.bcf.gz \
 --threads $numthreads --allow-overlaps -Oz \
 -o $outdir/filtered_consensus_ready_no-indels.vcf.gz

#you could uncomment this to remove the tmp files, but I prefer to remove by hand, in case something goes wrong with the script.
#rm $outdir/tmp-invariants.filtered.bcf.gz*
#rm $outdir/tmp-variants.filtered.bcf.gz*


#index final merged file
tabix $outdir/filtered_consensus_ready_no-indels.vcf.gz

