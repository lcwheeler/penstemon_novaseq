#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=astral_tree


cd $SLURM_SUBMIT_DIR

#path to ASTRAL
astral="/work/bs66/software/ASTRAL/Astral/astral.5.7.8.jar"

#path to window trees, catted into combined treefile
treefile="/work/bs66/dasanthera_novaseq/analysis/astral_trees/combined_windowtrees.tre"

#path to output file
outpath="/work/bs66/dasanthera_novaseq/analysis/astral_trees"


#run astral. Using single ML estimates of gene trees rather than bootstrap resampling.
#wants unrooted gene trees (which ours are -- the "rooting" in IQtree is just for visualization, and they are not rooted in the treefile)
#the branch annotations are for full annotations
java -jar $astral \
 --input $treefile 
 --outgroup mon_61-7_S440 \
 --branch-annotate 2 \
 --output $outpath/astral_out.tre

