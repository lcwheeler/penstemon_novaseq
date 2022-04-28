#!/bin/sh

#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p wessinger-48core
#SBATCH --job-name=testrun_fastp

cd $SLURM_SUBMIT_DIR


#path to fastp
fastpdir="/work/bs66/software"

#path to illumina reads
reads="/work/bs66/dasanthera_novaseq/merged_reads"


#fastp for loop
for r1in in $reads/*_R1_001.fastq.gz; 
do
    r2in="${r1in/R1_001.fastq.gz/R2_001.fastq.gz}"
    r1out="${r1in##*/}"
    r2out="${r1out/R1_001.fastq.gz/R2_001.fastq.gz}"
    $fastpdir/./fastp -i "$r1in" -I "$r2in" -m --merged_out "${r1out/merged_L001_R1_001.fastq.gz/trimmed.fastq.gz}" --include_unmerged -x -c -w 16 -h "${r1out/merged_L001_R1_001.fastq.gz/html}" -j "${r1out/merged_L001_R1_001.fastq.gz/json}"
done

