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

#path to directory holding subdirectories -- one for each scaffold
#also set up variable to infer trees in each scaffold directory
aligndir="/work/bs66/dasanthera_novaseq/analysis/genetree_infiles"


indirs=($aligndir/*/)
scafdir="${indirs[$SLURM_ARRAY_TASK_ID]}"
outdir=${scafdir/genetree_infiles/genetree_outfiles}



#submit iqtree job for each scaffold
iqtree -s $scafdir \
 --seqtype DNA \
 -o bws_mon_61-7_S440 \
 --verbose #\
# --prefix $outdir



