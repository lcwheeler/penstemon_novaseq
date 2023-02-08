#!/bin/sh

#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p defq
#SBATCH --job-name=aBSREL
#SBATCH --output=slurm-run-aBSREL_%j.out

# This script uses hyphy to fit the aBSREL model to codon alignment data for *all* branches 
# We use the species tree for the analyses for each gene


cd $SLURM_SUBMIT_DIR


#path to CDS.bed file from the annotation, and desired outfile
alndir="/work/lw74/habro/selection_analyses/CDS_reduced_BED_pal2nal_alignments"
mkdir -p "/work/lw74/habro/selection_analyses/CDS_reduced_BED_aBSREL_outputs_AllBranches"
outdir="/work/lw74/habro/selection_analyses/CDS_reduced_BED_aBSREL_outputs_AllBranches"
absrel_tree="/work/lw74/habro/selection_analyses/Habroanthus_WithOutgroups_ReducedBED_CDS_astral_out.reformat.tre"


# iterate over the codon aligments and run the aBSREL model in hyphy with species tree
# loop over the 3rd quarter of the alignments 
flist=$(ls $alndir/*.pal2nal | tail --lines=12191 | head --lines=6095)

for fas in $flist
do
	genename=$(basename $fas)
	genename=$(echo $genename | sed 's/.pal2nal//g')
	hyphy absrel --alignment $fas --tree $absrel_tree --CPU 16 --output $outdir/$genename".AllBranches.absrel.JSON" > $outdir/$genename".AllBranches.absrel.txt"
done

