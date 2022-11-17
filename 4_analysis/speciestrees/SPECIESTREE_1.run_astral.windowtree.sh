#!/bin/sh

#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p defq
#SBATCH --job-name=astral_run
#SBATCH --output=astral_run%j.out


cd $SLURM_SUBMIT_DIR


# Path to astral tree outfile

#mkdir -p "/work/lw74/habro/consensus_alignments/CDS_fastas/species_trees_habro"
outdir="/work/lw74/habro/species_trees_habro"
outfile="/work/lw74/habro/species_trees_habro/Habroanthus_windows_WithOutgroups_astral_out.tre"


# Path to file with concatenated gene trees for the genomic windows
trefile="/work/lw74/habro/consensus_alignments_WithOutgroups/individual_fullgenome_fastas/genetree_infiles/full_set_gene_trees/tidy_trees/Habroanthus_window_trees_concatenated.trees"


# Use Astral to generate a species tree using the full CDS gene trees
 
java -jar /work/lw74/software_installs/ASTRAL-5.7.1/Astral/astral.5.7.1.jar --branch-annotate 2 -i $trefile -o $outfile 2>$outdir/"habro_astral_windows_WithOutroupgs_out.log"

