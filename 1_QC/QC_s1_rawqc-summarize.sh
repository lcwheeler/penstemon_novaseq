#!/bin/sh

#SBATCH -N 1
#SBATCH -n 20
#SBATCH -p wessinger-48core
#SBATCH --job-name=QC_s1

cd $SLURM_SUBMIT_DIR


#source appropriate environments to enable use of conda installs through batch submission
source /home/lw74/.bashrc
source /home/lw74/.bash_profile

#activate conda environment with QC packages installed
#conda activate QC


################
#fastqc+multiqc#
################

rawreads="/work/lw74/habro"

fastqc_outdir_rawreads="/work/lw74/habro/fastqc_rawreads"
rm -r $fastqc_outdir_rawreads
mkdir $fastqc_outdir_rawreads

#move to directory with reads (raw reads directory)
cd $rawreads

#store fastq files as array
files=(*.fastq.gz)

#perform fastqc -- distinction from for loop is this can process -t files simultaneously
fastqc "${files[@]}" -t 19 -o $fastqc_outdir_rawreads

#summarize results with multiqc
multiqc $fastqc_outdir_rawreads -o $fastqc_outdir_rawreads
