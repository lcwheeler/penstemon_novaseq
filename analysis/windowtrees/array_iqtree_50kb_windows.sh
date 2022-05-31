#!/bin/sh

#SBATCH -N 1
#SBATCH -n 10
#SBATCH -p wessinger-48core
#SBATCH --job-name=test_iqtree

##NOTE: this is an arrayed batch script. It requires special syntax to submit jobs.
#I have 12 total scaffolds, and because counting starts at 0:
#submit with sbatch --array [0-11] array_mapping_pipe.sh


cd $SLURM_SUBMIT_DIR

#load conda env. Need numpy, matplotlib, biopython, (others?). installed.
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc

#paths to important files and directories
aligndir="/work/bs66/dasanthera_novaseq/50kb_alignments/"

#identify input files
cd $aligndir
base_list=(*_50kb_window)
current_dir="${base_list[$SLURM_ARRAY_TASK_ID]}"
cd $current_dir

for i in *.fa;
do
    iqtree -s $i -T 10 -o bws_mon_61-7_S440_mapped_marked_filtered.bam
done


