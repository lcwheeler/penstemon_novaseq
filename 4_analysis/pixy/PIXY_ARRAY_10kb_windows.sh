#!/bin/sh

#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p wessinger-48core
#SBATCH --job-name=pixy_run

# This runs pixy on the filtered all-states VCF file to calculate dXY statistics
# This is an array script that must be submitted using sbatch --array [0-n], where n = # linkage groups - 1
# The array setup allows each chromosome alignment to be run independently in parallel


# INPUTS: inputs to this script are the filtered full-genome VCF and pixy pop file


cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

# Source appropriate environments to enable use of conda installs through batch submission
source /home/lw74/.bashrc
source /home/lw74/.bash_profile

# Activate the pixy conda environment
conda activate pixy_env

# Input files
invcf="/work/lw74/habro/VCFs_WithOutgroups/filtered_consensus_ready_no-indels.vcf.gz"
popfile="/work/lw74/habro/pixy_analyses/popfile_pixy.txt"


# Make array for specifying the linkage groups
faidxfile="/work/lw74/refGenome/Pbar.2022.LG.fa.fai"
chromlist=$(awk '{print $1}' $faidxfile)
chromlist=(${chromlist[@]})
chromname="${chromlist[$SLURM_ARRAY_TASK_ID]}"


# Set the output directory, specificy sub-dir name to match chrom name
mkdir -p "/work/lw74/habro/pixy_analyses/pixy_outputs"
outdir="/work/lw74/habro/pixy_analyses/pixy_outputs/$chromname"


##########
#ANALYSIS#
##########

# Run pixy on the full Habroanthus dataset

pixy --stats pi fst dxy \
 --vcf $invcf \
 --populations $popfile \
 --chromosomes $chromname \
 --window_size 10000 \
 --n_cores 4 \
 --output_folder $outdir \
 --output_prefix $chromname
