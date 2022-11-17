#!/bin/sh

#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p wessinger-48core
#SBATCH --job-name=combine_new_seqs
#SBATCH --output=slurm-combine_new_seqs%j.out


cd $SLURM_SUBMIT_DIR


#path to output directory for the gene trees

mkdir -p "/work/lw74/habro/jonesii_side_analysis/consensus_alignments/CDS_fastas/combined_full_CDS"
outdir="/work/lw74/habro/jonesii_side_analysis/consensus_alignments/CDS_fastas/combined_full_CDS"


#path to directory with alignments (fasta format) for each CDS
alndir="/work/lw74/habro/jonesii_side_analysis/consensus_alignments/CDS_fastas/matched_subset_seqs"

#path to the directory where the subset of rearranged, filtered jonesii CDS files are
indir="/work/lw74/habro/jonesii_side_analysis/consensus_alignments/CDS_fastas/subset_seqs"



# Loop over all the CDS alignments and run FastTree

flist=$(cat $SLURM_SUBMIT_DIR/matched_filenames.txt)

for infile in $flist
do
	
    cat $alndir/$infile >> $outdir/$infile
    cat $indir/$infile >> $outdir/$infile

done

