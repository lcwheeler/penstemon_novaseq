#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=filter_passed_relocate

cd $SLURM_SUBMIT_DIR

indir="/work/lw74/habro/jonesii_side_analysis/consensus_alignments/CDS_fastas/filtered_consensus_fasta"
outdir="/work/lw74/habro/jonesii_side_analysis/consensus_alignments/CDS_fastas/passed_filter"
mkdir -p "/work/lw74/habro/jonesii_side_analysis/consensus_alignments/CDS_fastas/passed_filter"

flist="/work/lw74/habro/jonesii_side_analysis/consensus_alignments/CDS_fastas/passed-files.txt"

cd $indir

mv -t $outdir $(< $flist)
