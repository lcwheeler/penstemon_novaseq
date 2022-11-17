#!/bin/sh

#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p wessinger-48core
#SBATCH --job-name=fasttree_genetrees

##NOTE: this is an arrayed batch script. The array is 0-n, where n is the number of scaffolds-1. For Habroanthus analyses, this is sbatch --array [0-10]


cd $SLURM_SUBMIT_DIR

#load conda env. Need iqtree installed.
source /home/lw74/.bashrc
source /home/lw74/.bash_profile


#path to directory holding subdirectories -- one for each scaffold
aligndir="/work/lw74/habro/consensus_alignments_WithOutgroups/individual_fullgenome_fastas/genetree_infiles"


#set up variable to infer trees in each scaffold directory
indirs=($aligndir/*/)
inscaf="${indirs[$SLURM_ARRAY_TASK_ID]}"


#change to the input directory
cd $inscaf

#for loop to go through each alignment, estimate gene tree using FastTree

for f in *.fa
do
        FastTree -gtr -nt $f > $f".tre"
done

#list outdir file and move files
outdir=${inscaf/genetree_infiles/genetree_outfiles}
mv outfile* $outdir
