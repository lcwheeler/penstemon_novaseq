#!/bin/sh

#SBATCH -N 1
#SBATCH -n 12
#SBATCH -p wessinger-48core
#SBATCH --job-name=astral_run
#SBATCH --output=astral_run%j.out


cd $SLURM_SUBMIT_DIR


# Path to astral tree outfile

mkdir -p "/work/lw74/habro/consensus_alignments/CDS_fastas/species_trees_habro"
outdir="/work/lw74/habro/consensus_alignments/CDS_fastas/species_trees_habro"
outfile="/work/lw74/habro/consensus_alignments/CDS_fastas/species_trees_habro/Habroanthus_astral_out.tre"


# Path to file with concatenated gene trees for the full CDS
trefile="/work/lw74/habro/consensus_alignments/CDS_fastas/gene_trees_full_CDS/Habroanthus_concatenated_gene_trees.tree"


# Use Astral to generate a species tree using the full CDS gene trees
 
java -jar /work/lw74/software_installs/ASTRAL-5.7.1/Astral/astral.5.7.1.jar -i $trefile -o $outfile 2>$outdir/"habro_astral_out.log"
