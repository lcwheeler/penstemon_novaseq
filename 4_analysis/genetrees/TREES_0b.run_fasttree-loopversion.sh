#!/bin/sh

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p wessinger-48core
#SBATCH --job-name=fasstree_gene_trees
#SBATCH --output=slurm-fasttree%j.out


cd $SLURM_SUBMIT_DIR


#path to output directory for the gene trees

mkdir -p "/work/lw74/habro/consensus_alignments_WithOutgroups/CDS_fastas/gene_trees_full_CDS"
outdir="/work/lw74/habro/consensus_alignments_WithOutgroups/CDS_fastas/gene_trees_full_CDS"
concattree="/work/lw74/habro/consensus_alignments_WithOutgroups/CDS_fastas/gene_trees_full_CDS/Habroanthus_concatenated_gene_trees.tree"

#path to directory with alignments (fasta format) for each CDS
alndir="/work/lw74/habro/consensus_alignments_WithOutgroups/CDS_fastas/filtered_consensus_fasta/passed_files"


# Loop over all the CDS alignments and run FastTree

flist=$(ls *.cln.fa | head --lines=11827)

for infile in $flist
do
	
    # run FastTree on the CDS alignments to generate gene trees
    FastTree -gtr -nt $infile > $infile".tre"

done

# Concatenate the gene trees for Astral input and relocate to the out directory
cat *.tre > $concattree
mv *.tre $outdir

