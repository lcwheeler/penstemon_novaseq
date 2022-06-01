#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=astral_v1


cd $SLURM_SUBMIT_DIR

#path to ASTRAL
astral="/work/bs66/software/ASTRAL/Astral/astral.5.7.8.jar"

#path to window trees, catted into combined treefile
treefile="/work/bs66/dasanthera_novaseq/50kb_alignments/combined_50kb_trees.tre"

#path to output file
outpath="/work/bs66/dasanthera_novaseq/analysis_astral"


#run astral
java -jar $astral -i $treefile -a $outpath/namemap_astral_v1.txt --outgroup P_montanus -o $outpath/astral_tree_50kb_v1.tre
