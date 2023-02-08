#!/bin/sh

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p wessinger-48core
#SBATCH --job-name=conflict


#wessinger-48core
cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

########################
#source appropriate environments to enable use of conda installs through batch submission
source /home/lw74/.bashrc
source /home/lw74/.bash_profile


#activate conda environment (python 2.7) to run scripts
conda activate py27

#Location of species tree and gene tree files
sptree="/work/lw74/habro/species_trees_habro/Habroanthus_CDS_WithOutgroups_astral_out.tre"
gtree_dir="/work/lw74/habro/consensus_alignments_WithOutgroups/CDS_fastas/gene_trees_full_CDS"

#Define outfile prefix for outputs
prefix="CDS_WithGroups"

#*.tre
########################

#####
#Run the python program to analyze gene tree conflict#
#####


python /work/lw74/software_installs/gene_family_conflicts/scripts/CAnDI.py --mode n --species_tree $sptree --gene_folder $gtree_dir --outfile_prefix $prefix

# don't use --output_subtrees in this analysis
#mkdir -p sub_trees_output
#mv *.sub.tre sub_trees_output
