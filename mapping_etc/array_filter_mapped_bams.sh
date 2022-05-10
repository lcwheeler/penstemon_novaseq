#!/bin/sh

#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p wessinger-48core
#SBATCH --job-name=filter_mapped_bams

##NOTE: this is an arrayed batch script. It requires special syntax to submit job.
#I have 18 total samples, and because counting starts at 0:
#submit with sbatch --array [0-17] array_mapping_pipe.sh

cd $SLURM_SUBMIT_DIR

source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc

#paths to important files and directories
refgenome="/work/bs66/project_compare_genomes/filtered_davidsonii_genome_1mb.fasta"
mappedreads="/work/bs66/dasanthera_novaseq/mapped_marked_bams"
outdir="/work/bs66/dasanthera_novaseq/mapped_marked_filtered_bams"


#identify reads
cd $mappedreads
readlist=(*.bam)
bamread="${readlist[$SLURM_ARRAY_TASK_ID]}"

#1. quality score filter (-q 20), remove marked duplicates (-F 0x400)

samtools view -@ 4 -h -q 20 -F 0x400 $bamread | samtools sort -@ 4 -o $outdir/"${bamread/mapped_marked.bam/mapped_marked_filtered.bam}"

#index reads
samtools index -b $outdir/"${bamread/mapped_marked.bam/mapped_marked_filtered.bam}"

