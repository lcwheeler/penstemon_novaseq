#!/bin/sh

#SBATCH -N 1
#SBATCH -n 20
#SBATCH -p wessinger-48core
#SBATCH --job-name=fastqc_trimmed

cd $SLURM_SUBMIT_DIR

#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

#activate conda environment with fastqc and multiqc installed
conda activate QC


#path to filtered illumina reads
reads="/work/bs66/dasanthera_novaseq/filtered_reads"
outdir="/work/bs66/dasanthera_novaseq/fastqc_filtered"

#move to directory with reads
cd $reads

#store fastq files as array
files=(*.fastq.gz)

#perform fastqc -- distinction from for loop is this can process -t files simultaneously
fastqc "${files[@]}" -t 20 -o $outdir


