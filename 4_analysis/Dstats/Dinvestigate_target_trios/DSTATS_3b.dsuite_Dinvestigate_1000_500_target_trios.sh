#!/bin/sh

#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p defq
#SBATCH --job-name=Dstat_dinv
#SBATCH --output=slurm-Dsuite_Dinvestigate_%j.out
#SBATCH --error=slurm-Dsuite_Dinvestigate_%j.error


cd $SLURM_SUBMIT_DIR

# Running this analysis with a set of target triplets that include only the monophyletic hummingbird and bees species from sister pairs
# e.g. barbatus, eatonii

#######
#SETUP#
#######

#source appropriate environments to enable use of conda installs through batch submission
source /home/lw74/.bashrc
source /home/lw74/.bash_profile


#set Dsuite variable to path, and for f-branch test, the path to directory with dtools.py
dtools="/work/lw74/software_installs/Dsuite/utils"


#specify invcfdir, invcf, outdir, faidx file, tree file
invcf="/work/lw74/habro/VCFs_WithOutgroups/filtered_consensus_ready_no-indels.vcf.gz"
faidxfile="/work/lw74/refGenome/Pbar.2022.LG.fa.fai"
test_trios="/work/lw74/habro/Dstats_analyses/test_trios_target_trios.txt"
popset="/work/lw74/habro/Dstats_analyses/popset_dtrios_target_trios.txt"
dstat_outdir="/work/lw74/habro/Dstats_analyses/dsuite_outputs/dinvestigate_results_target_trios"

mkdir -p "/work/lw74/habro/Dstats_analyses/dsuite_outputs/dinvestigate_results_target_trios"


##########
#ANALYSIS#
##########


#perform investigate analyses -- try large windows first, and on whole genome 
Dsuite Dinvestigate -w 1000,500 $invcf $popset $test_trios

mv *_1000_500.txt $dstat_outdir

