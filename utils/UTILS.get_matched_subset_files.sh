#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=filter_passed_relocate

cd $SLURM_SUBMIT_DIR

indir="/work/lw74/habro/consensus_alignments_WithOutgroups/CDS_fastas/filtered_consensus_fasta/passed_files/"
reffile="/work/lw74/habro/jonesii_side_analysis/consensus_alignments/CDS_fastas/subset_filenames.txt"

outdir="/work/lw74/habro/jonesii_side_analysis/consensus_alignments/CDS_fastas/matched_subset_seqs"
mkdir -p "/work/lw74/habro/jonesii_side_analysis/consensus_alignments/CDS_fastas/matched_subset_seqs"


cd $indir

xargs -a $reffile cp -t $outdir
