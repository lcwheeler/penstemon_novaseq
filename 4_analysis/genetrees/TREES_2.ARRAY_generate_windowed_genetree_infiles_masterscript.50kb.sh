#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=generate_genetree_infiles

#this is a masterscript to generate sliding window trees 
#This is an array script. Array is 0 to n-1, where n is number of scaffolds in scaflist


cd $SLURM_SUBMIT_DIR


#source bash profile and activate conda environment with python3 as base python env
source /home/lw74/.bashrc
source /home/lw74/.bash_profile


#paths to important scripts and directories
pythonscript="/work/lw74/habro/consensus_alignments_WithOutgroups/individual_fullgenome_fastas/TREES.create_fasta_window_alignments-MOD.py"
scaffold_dir="/work/lw74/habro/consensus_alignments_WithOutgroups/individual_fullgenome_fastas/rearranged_consensus_seqs_fullgenomes"
outdir="/work/lw74/habro/consensus_alignments_WithOutgroups/individual_fullgenome_fastas/50kb_windows"


#list of scaffold names
scaflist=("LG01" "LG02" "LG03" "LG04" "LG05" "LG06" "LG07" "LG08" "LG09" "LG10" "LG11")
inscaf="${scaflist[$SLURM_ARRAY_TASK_ID]}"

#run for loop for python script to generate gene tree files for each scaffold
#syntax: python3 scriptname input_fasta windowsize prefix

mkdir -p $outdir
python $pythonscript $scaffold_dir/$inscaf.rearranged.fa 50000 0.75 $outdir/$inscaf "_consensus_fullgenome_" 

