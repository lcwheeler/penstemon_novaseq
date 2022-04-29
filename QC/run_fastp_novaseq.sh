#!/bin/sh

#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p wessinger-48core
#SBATCH --job-name=testrun_fastp

cd $SLURM_SUBMIT_DIR


#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

#activate conda environment with fastqc and multiqc installed
conda activate QC

#path to illumina reads
reads="/work/bs66/dasanthera_novaseq/merged_reads"


#fastp for loop
for r1in in $reads/*_R1_001.fastq.gz; 
do
    r2in="${r1in/R1_001.fastq.gz/R2_001.fastq.gz}"
    r1out="${r1in##*/}"
    r2out="${r1out/R1_001.fastq.gz/R2_001.fastq.gz}"
    fastp -i "$r1in" -I "$r2in" --out1 "${r1out/merged_L001_R1_001.fastq.gz/trimmed_L001_R1_001.fastq.gz}" --out2 "${r1out/merged_L001_R1_001.fastq.gz/trimmed_L001_R2_001.fastq.gz}" --unpaired1 "${r1out/merged_L001_R1_001.fastq.gz/unpaired_L001_R1_001.fastq.gz}" --unpaired2 "${r1out/merged_L001_R1_001.fastq.gz/unpaired_L001_R2_001.fastq.gz}" -x -c -w 16 -h "${r1out/merged_L001_R1_001.fastq.gz/html}" -j "${r1out/merged_L001_R1_001.fastq.gz/json}"
done

