#!/bin/sh

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p wessinger-48core
#SBATCH --job-name=filter_rearranged_consensus

#wessinger-48core
cd $SLURM_SUBMIT_DIR


# Run the custom Python script to filter fasta files by percentage missing data

python missing-data-filter-args.py -per_n 10 -fpattern "*.rearranged.fa" -rm_ext ".rearranged.fa" -save_ext ".rearranged.cln.fa" -num_samples 42 -per_m 100 -id_split_string "full_"


# Relocate the filtered output files to a new directory

mkdir -p filtered_consensus_fasta
mv *.cln.fa filtered_consensus_fasta
mv missing-data-counts.csv filtered_consensus_fasta
mv removed-seqs-stats.csv filtered_consensus_fasta
mv total-sequence-length-data.csv filtered_consensus_fasta
mv passed-files.txt filtered_consensus_fasta
mv filtered_consensus_fasta ../
