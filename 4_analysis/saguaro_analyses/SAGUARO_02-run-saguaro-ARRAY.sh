#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=saguaro_run

# This runs Saguaro on the dataset after conversion of multi-fasta alignment to binary features
# This is an array script that must be submitted using sbatch --array [0-n], where n = # input files - 1
# The array setup allows each chromosome alignment to be run independently in parallel


# INPUTS: inputs to this script are the binary Saguaro feature files that were converted from multi-fasta format


cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

#source appropriate environments to enable use of conda installs through batch submission
source /home/lw74/.bashrc
source /home/lw74/.bash_profile


# Make array for specifying the linkage groups
featuredir="/work/lw74/habro/saguaro_analyses"
chromlist=$(ls $featuredir/*.features)
chromlist=(${chromlist[@]})
chromfile="${chromlist[$SLURM_ARRAY_TASK_ID]}"
chromname="$(basename $chromfile | sed 's/.features//g')"


# Set the output directory, specificy sub-dir name to match chrom name
mkdir -p "/work/lw74/habro/saguaro_analyses/saguaro_outputs"
outdir="/work/lw74/habro/saguaro_analyses/saguaro_outputs/$chromname"


##########
#ANALYSIS#
##########

# Run Saguaro on Habroanthus data 

Saguaro -f $chromfile -o $outdir -iter 40
