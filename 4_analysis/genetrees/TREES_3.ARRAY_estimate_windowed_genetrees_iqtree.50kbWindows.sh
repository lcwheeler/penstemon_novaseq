#!/bin/sh

#SBATCH -N 1
#SBATCH -n 6
#SBATCH -p defq
#SBATCH --job-name=iqtree_genetrees

##NOTE: this is an arrayed batch script. The array is 0-n, where n is the number of scaffolds-1. For P. barbatus genome, this is sbatch --array [0-8]


cd $SLURM_SUBMIT_DIR

#load conda env. Need iqtree installed.
source /home/lw74/.bashrc
source /home/lw74/.bash_profile


#path to directory holding subdirectories -- one for each scaffold
alndir="/work/lw74/habro/consensus_alignments_WithOutgroups/individual_fullgenome_fastas/50kb_windows"

mkdir -p "/work/lw74/habro/consensus_alignments_WithOutgroups/individual_fullgenome_fastas/IQtree_50kb_window_trees"
outdir="/work/lw74/habro/consensus_alignments_WithOutgroups/individual_fullgenome_fastas/IQtree_50kb_window_trees"
concat_tree="/work/lw74/habro/consensus_alignments_WithOutgroups/individual_fullgenome_fastas/IQtree_50kb_window_concatenated_genetrees.tre"
concat_tree_log="/work/lw74/habro/consensus_alignments_WithOutgroups/individual_fullgenome_fastas/IQtree_50kb_window_concatenated_genetrees_order.txt"


#set up variable to infer trees in each scaffold directory
#indirs=($alndir/LG*)
scaflist=(LG01 LG02 LG03 LG04 LG05 LG06 LG07 LG08 LG09 LG10 LG11)
scaflist=$(scaflist[@])
inscaf="${scaflist[$SLURM_ARRAY_TASK_ID]}"


#change to the input directory
cd $alndir

#for loop to go through each alignment, estimate subs. model, and est. tree
#Maybe there is better way but I can't figure out other than for loop
#the ARRAY number changes the scaffold for which this occurs
#"rooting" at P. palmeri (CAW98_1_S441), but this is still an unrooted tree

for i in $inscaf*.fa
do
	iqtree -s $i \
	 --seqtype DNA \
	 -m MFP \
	 -o CAW98_1_S441 \
	 -T 6 \
	 --prefix outfile_$i
done

#list outdir file and move files
#mv outfile* $outdir

# concatenate trees
#cd $outdir
#ls *.treefile > $concat_tree_log
#cat *.treefile > $concat_tree
