#!/bin/sh

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p wessinger-48core
#SBATCH --job-name=filter_rearranged_consensus

#wessinger-48core
cd $SLURM_SUBMIT_DIR

indir="/work/lw74/habro/consensus_alignments_WithOutgroups_ReducedBED/CDS_fastas/rearranged_consensus_seqs_CDS_full"
outdir="/work/lw74/habro/consensus_alignments_WithOutgroups_ReducedBED/CDS_fastas/filtered_consensus_fasta"
passdir="/work/lw74/habro/consensus_alignments_WithOutgroups_ReducedBED/CDS_fastas/passed_filtered_consensus"

# Path to the python script
script="/work/lw74/habro/consensus_alignments_WithOutgroups_ReducedBED/CDS_fastas/missing-data-filter-args.py"

cd $indir

date

# Run the custom Python script to filter fasta files by percentage missing data

python $script -per_n 25 -fpattern "*.rearranged.fa" -rm_ext ".rearranged.fa" -save_ext ".rearranged.cln.fa" -num_samples 45 -id_split_string "full_"

# Relocate the filtered output files to a new directory

mkdir -p $outdir
mkdir -p passdir
mv *.cln.fa $outdir
mv missing-data-counts.csv $outdir
mv removed-seqs-stats.csv $outdir
mv total-sequence-length-data.csv $outdir
mv passed-files.txt $outdir

# Move to the filtered dir and relocate only the passed files to a new dir
cd $outdir
mv -t $passdir $(< passed-files.txt)

date
