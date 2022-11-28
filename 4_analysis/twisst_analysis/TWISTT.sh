#!/bin/sh

#SBATCH -N 1
#SBATCH -n 12
#SBATCH -p defq
#SBATCH --job-name=twisst_run


cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

#source appropriate environments to enable use of conda installs through batch submission
source /home/lw74/.bashrc
source /home/lw74/.bash_profile


#specify invcfdir, invcf, outdir, faidx file, tree file
twisst_loc="/work/lw74/software_installs/twisst"
input_trees="/work/lw74/habro/twisst_analysis/genomic_windows_analysis/Habroanthus_window_trees_concatenated.trees"
outfile="/work/lw74/habro/twisst_analysis/genomic_windows_analysis/output_windows_astral_tree.weights.csv.gz"
groupsfile="/work/lw74/habro/twisst_analysis/genomic_windows_analysis/groups.tsv"


##########
#ANALYSIS#
##########

# Use twisst to calculate the topology weights for the Habroanthus gene tree data
# Group memberships are specified in the python script call e.g. -g barbatus 

#python $twisst_loc/twisst.py -t $input_trees -w $outfile -g specios -g laevis -g eatonii -g pseudoputus -g putus -g labrosus -g comarrhenus -g virgatus_virgatus -g neomexicanus -g barbatus -g cardinalis -g glaber_alpinus --method complete --groupsFile $groupsfile
python $twisst_loc/twisst.py -t $input_trees -w $outfile -g specios -g laevis -g eatonii -g pseudoputus -g labrosus -g virgatus_virgatus -g barbatus -g cardinalis --method complete --groupsFile $groupsfile


