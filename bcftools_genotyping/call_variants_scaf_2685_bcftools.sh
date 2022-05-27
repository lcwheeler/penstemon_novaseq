#!/bin/sh

#SBATCH -N 1
#SBATCH -n 10
#SBATCH -p wessinger-48core
#SBATCH --job-name=variants_bcftools


cd $SLURM_SUBMIT_DIR

source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc

#paths to important files and directories
refgenome="/work/bs66/project_compare_genomes/annot_Pdavidsonii_genome.1mb.fasta"
mapped_filtered_reads="/work/bs66/dasanthera_novaseq/mapped_marked_filtered_bams"
outdir="/work/bs66/dasanthera_novaseq/called_variants_bcftools"
scaffold="scaffold_2685"


#mpileup
#-Ou: output, (u) uncompressed added because piping
#-a: annotate (FORMAT/AD) allelic depth

#call
#-v: variants only, multiallelic caller, 
#-m: multiallelic caller
#-O: output format (v) uncompressed vcf

#2. pileup
#3. call variants

bcftools mpileup --threads 10 -Ou -I -a FORMAT/AD --max-depth 100 -r $scaffold -f $refgenome $mapped_filtered_reads/*.bam | bcftools call --threads 10 -m -Ou | bcftools filter -g 5 -S . -o $outdir/alignment.$scaffold.vcf
