#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=saguaro_convert

# This script converts alignment data for the Habroanthus species into Saguaruo binary format
# This is an array script that must be submitted using sbatch --array [0-n], where n = # input files - 1
# The array setup allows each chromosome alignment to be run independently in parallel

# INPUTS: input to the script are the rearranged full-genome fasta files extracted from the full filtered
# VCF. These should have one entry per sample and be separated by chromsome (e.g. LG01.rearranged.fa)

cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

#source appropriate environments to enable use of conda installs through batch submission
source /home/lw74/.bashrc
source /home/lw74/.bash_profile


# Make array for specifying the linkage groups 
mafdir="/work/lw74/habro/consensus_alignments_WithOutgroups/individual_fullgenome_fastas/rearranged_consensus_seqs_fullgenomes"
chromlist=$(ls $mafdir/*.rearranged.fa)
chromlist=(${chromlist[@]})
chromfile="${chromlist[$SLURM_ARRAY_TASK_ID]}"
chromname="$(basename $chromfile | sed 's/.rearranged.fa//g')"


# Set the input files
namesfile="/work/lw74/habro/saguaro_analyses/names.txt"


# Set the output files
outdir="/work/lw74/habro/saguaro_analyses"
outbin="$outdir/$chromname.saguaro.features"


##########
#ANALYSIS#
##########

# Run Saguaro on Habroanthus data with subset of species

Fasta2HMMFeature -i $chromfile -o $outbin -nosame
