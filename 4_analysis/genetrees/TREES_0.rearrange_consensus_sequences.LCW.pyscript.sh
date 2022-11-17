#!/bin/sh

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p wessinger-48core
#SBATCH --job-name=rearrange_consensus

#wessinger-48core
cd $SLURM_SUBMIT_DIR

# Run custom python script to rearrange consensus sequences
python rearrange-consensus-seqs.py

mkdir -p rearranged_consensus_seqs_CDS_full
mv *.rearranged.fa rearranged_consensus_seqs_CDS_full
