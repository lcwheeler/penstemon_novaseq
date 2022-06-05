#!/bin/sh

#SBATCH -N 1
#SBATCH -n 5
#SBATCH -p wessinger-48core
#SBATCH --job-name=cf_50kb_nomap


cd $SLURM_SUBMIT_DIR

#load conda env with iqtree installed
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc

#paths to important files and directories
reftree="/work/bs66/dasanthera_novaseq/analysis_astral/astral_tree_50kb_v1_no-namemap_rooted.tre"
sourcetrees="/work/bs66/dasanthera_novaseq/50kb_alignments/combined_50kb_trees.tre"
outdir="/work/bs66/dasanthera_novaseq/analysis_concordance_factors_50kb"


cd $outdir
iqtree -t $reftree --gcf $sourcetrees --cf-verbose --df-tree --prefix no_namemap -T 5

