#!/bin/sh

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p wessinger-48core
#SBATCH --job-name=twisst_run


cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

#source appropriate environments to enable use of conda installs through batch submission
source /home/lw74/.bashrc
source /home/lw74/.bash_profile


#specify twisst input and output files
twisst_loc="/work/lw74/software_installs/twisst"
input_trees="/work/lw74/habro/twisst_analysis/CDS_target_trios/IQtree_CDS_concatenated_genetrees.tre"
groupsfile="/work/lw74/habro/twisst_analysis/CDS_target_trios/groups.tsv"
outdir="/work/lw74/habro/twisst_analysis/CDS_target_trios"


##########
#ANALYSIS#
##########

# Use twisst to calculate the topology weights for the Habroanthus gene tree data
# Groups included in the analysis are specified in the python script call e.g. -g barbatus 
# Individual samples are assigned to groups in the groups.tsv file
# Run a set of trio comparisons matched to those that I ran the fDM stats on with Dsuite
# Use P. palmeri as the outgroup for all comparisons

# Run Twisst for tests of introgression between barbatus and eatonii

python $twisst_loc/twisst.py -t $input_trees \
 -w $outdir/virgatus_barbatus_eatonii.iqree.weights.txt.gz \
 -g virgatus -g barbatus -g eatonii -g palmeri --outgroup palmeri \
 --method complete --groupsFile $groupsfile

python $twisst_loc/twisst.py -t $input_trees \
 -w $outdir/laevis_eatonii_barbatus.iqtree.weights.txt.gz \
 -g laevis -g eatonii -g barbatus -g palmeri --outgroup palmeri \
 --method complete --groupsFile $groupsfile
 
 
# Run Twisst for tests of introgression between barbatus and labrosus

python $twisst_loc/twisst.py -t $input_trees \
 -w $outdir/virgatus_barbatus_labrosus.iqtree.weights.txt.gz \
 -g virgatus -g barbatus -g labrosus -g palmeri --outgroup palmeri \
 --method complete --groupsFile $groupsfile

python $twisst_loc/twisst.py -t $input_trees \
 -w $outdir/putus_labrosus_barbatus.iqtree.weights.txt.gz \
 -g putus -g labrosus -g barbatus -g palmeri --outgroup palmeri \
 --method complete --groupsFile $groupsfile



# Run Twisst for tests of introgression between barbatus and cardinalis

python $twisst_loc/twisst.py -t $input_trees \
 -w $outdir/virgatus_barbatus_cardinalis.iqtree.weights.txt.gz \
 -g virgatus -g barbatus -g cardinalis -g palmeri --outgroup palmeri \
 --method complete --groupsFile $groupsfile



# Run Twisst for tests of introgression between eatonii and labrosus

python $twisst_loc/twisst.py -t $input_trees \
 -w $outdir/laevis_eatonii_labrosus.iqtree.weights.txt.gz \
 -g laevis -g eatonii -g labrosus -g palmeri --outgroup palmeri \
 --method complete --groupsFile $groupsfile

python $twisst_loc/twisst.py -t $input_trees \
 -w $outdir/putus_labrosus_eatonii.iqtree.weights.txt.gz \
 -g putus -g labrosus -g eatonii -g palmeri --outgroup palmeri \
 --method complete --groupsFile $groupsfile



# Run Twisst for tests of introgression between eatonii and cardinalis

python $twisst_loc/twisst.py -t $input_trees \
 -w $outdir/laevis_eatonii_cardinalis.iqtree.weights.txt.gz \
 -g laevis -g eatonii -g cardinalis -g palmeri --outgroup palmeri \
 --method complete --groupsFile $groupsfile


# Run Twisst for tests of introgression between cardinalis and labrosus

python $twisst_loc/twisst.py -t $input_trees \
 -w $outdir/putus_labrosus_cardinalis.iqtree.weights.txt.gz \
 -g putus -g labrosus -g cardinalis -g palmeri --outgroup palmeri \
 --method complete --groupsFile $groupsfile












