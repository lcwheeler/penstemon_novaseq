#!/bin/sh

#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p wessinger-48core
#SBATCH --job-name=concordance


#wessinger-48core
cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

########################
#source appropriate environments to enable use of conda installs through batch submission
source /home/lw74/.bashrc
source /home/lw74/.bash_profile


#Location of species tree and gene tree files
species_treefile="/work/lw74/habro/species_trees_habro/Habroanthus_windows_WithOutgroups_astral_out.tre"
loci_treefile="/work/lw74/habro/consensus_alignments_WithOutgroups/individual_fullgenome_fastas/genetree_infiles/full_set_gene_trees/tidy_trees/Habroanthus_window_trees_concatenated.trees"


#Define outfile prefix for outputs
prefix="genomic_windows_WithOutgroups_concord"

########################

#####
#Run IQtree to calculate concordance factors
#####


iqtree2 -t $species_treefile --gcf $loci_treefile --prefix $prefix
