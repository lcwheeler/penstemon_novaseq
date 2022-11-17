#!/bin/sh

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p wessinger-48core
#SBATCH --job-name=makebed

###NOTE: this is set up as a batch array script. The array number is the number of values of K, -1

### This script converts a VCF input file to a plink bed file and then filters SNPs based on distance

cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

#source appropriate environments to enable use of conda installs through batch submission
source /home/lw74/.bashrc
source /home/lw74/.bash_profile

#vcfdir="/work/lw74/habro/VCFs_WithOutgroups"
invcf="/work/lw74/habro/VCFs_WithOutgroups/filtered_consensus_ready_no-indels.vcf.gz"
fbn=$(basename $invcf .vcf.gz)

# Convert the whole-genome filtered VCF to a plink bed file
#plink --vcf $invcf --make-bed --out $fbn --double-id --allow-extra-chr


# Thin the plink bed file by n bp separation between sites
bedfile=$fbn".bed"
bimfile=$fbn".bim"
famfile=$fbn".fam"

plink --bed $bedfile --bim $bimfile --fam $famfile --bp-space 1000 --make-bed --out $fbn"_thinned_1000bp" --double-id --allow-extra-chr


# Reformat chromosome names to integers in the bim file to make compatible with Admixture

cat $fbn"_thinned_1000bp".bim | sed 's/LG0//g' | sed 's/LG//g' > $fbn"_thinned_1000bp".bim.rfmt
mv $fbn"_thinned_1000bp".bim.rfmt > $fbn"_thinned_1000bp".bim
