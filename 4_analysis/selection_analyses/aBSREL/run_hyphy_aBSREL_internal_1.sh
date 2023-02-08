#!/bin/sh

#SBATCH -N 1
#SBATCH -n 10
#SBATCH -p wessinger-48core
#SBATCH --job-name=aBSREL
#SBATCH --output=slurm-generate-aBSREL.out

# This script uses hyphy to fit the aBSREL model to codon alignment data for all internal branches 
# We use the species tree for the analyses for each gene


cd $SLURM_SUBMIT_DIR


#path to CDS.bed file from the annotation, and desired outfile
alndir="/work/lw74/habro/selection_analyses/CDS_reduced_BED_pal2nal_alignments"
mkdir -p "/work/lw74/habro/selection_analyses/CDS_reduced_BED_aBSREL_outputs"
outdir="/work/lw74/habro/selection_analyses/CDS_reduced_BED_aBSREL_outputs"
absrel_tree="/work/lw74/habro/selection_analyses/Habroanthus_WithOutgroups_ReducedBED_CDS_astral_out.reformat.tre"


# iterate over the codon aligments and run the aBSREL model in hyphy with species tree

flist=$(ls $alndir/*.pal2nal | tail --lines=12191)

for fas in $flist
do
	genename=$(basename $fas)
	genename=$(echo $genename | sed 's/.pal2nal//g')
	hyphy absrel --alignment $fas --tree $absrel_tree --branches Internal --CPU 10 --output $outdir/$genename".absrel.JSON" > $outdir/$genename".absrel.txt"
done

