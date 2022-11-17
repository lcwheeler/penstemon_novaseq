#!/bin/sh

#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p wessinger-48core
#SBATCH --job-name=admixture_run
#SBATCH --output=slurm-admixture_%j.out
#SBATCH --error=slurm-admixture_%j.error


###NOTE: this is set up as a batch array script. The array number is the number of values of K, -1


cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

#source appropriate environments to enable use of conda installs through batch submission
source /home/lw74/.bashrc
source /home/lw74/.bash_profile

threads=4

# input files

bedfile="/work/lw74/habro/admixture_analysis/filtered_consensus_ready_no-indels_thinned_1000bp.bed"


##########
#ANALYSIS#
##########

# Create a list for the array input with different values of k
klist=(1 5 10 15 20 25 30 35 40)
klist=(${klist[@]})
K="${klist[$SLURM_ARRAY_TASK_ID]}"

# Loop version
#for K in 1 2 3 4 5; do admixture --cv hapmap3.bed $K | tee log${K}.out; done

# Run admixture on the bed input file with arrayed valued of K
admixture --cv $bedfile $K -j$threads | tee log${K}.out
