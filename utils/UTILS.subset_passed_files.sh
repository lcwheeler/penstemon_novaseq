#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=subset_passed_files

cd $SLURM_SUBMIT_DIR

indir="/work/lw74/habro/jonesii_side_analysis/consensus_alignments/CDS_fastas/passed_filter"
outdir="/work/lw74/habro/jonesii_side_analysis/consensus_alignments/CDS_fastas/subset_seqs"
mkdir -p "/work/lw74/habro/jonesii_side_analysis/consensus_alignments/CDS_fastas/subset_seqs"

cd $indir

ls *.cln.fa | shuf -n 2000 | xargs cp -t $outdir

cd $outdir
ls *.cln.fa > subset_filenames.txt
mv subset_filenames.txt $SLURM_SUBMIT_DIR
