#!/bin/sh

#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p wessinger-48core
#SBATCH --job-name=variants_bcftools

###NOTE: this is set up as a batch array script. Read carefully to understand how to submit properly.


cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

########################
#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

#activate conda environment with packages installed
#needs samtools (v1.15.1 works) and bcftools (v1.15.1 works)
conda activate mapping_etc

#load vcftools (v0.1.17 on cluster works)
#module load vcftools


#hard path to the reference genome and annotations file of interest
#mapped filtered reads should (1) not have a "/" at the end, and (2) be the same directory as output from script: mapping_pipeline_array_s1a-f.sh
refgenome="/work/bs66/project_compare_genomes/annot_Pdavidsonii_genome.1mb.fasta"
annotations_file="/work/bs66/project_compare_genomes/annot_Pdavidsonii_gffread.genes.bed"
mapped_filtered_reads="/work/bs66/dasanthera_novaseq/mapped_filtered_bams"


#tmp -- location of unfiltered input vcf
invcf="/work/bs66/dasanthera_novaseq/called_variants/unfiltered_gvcf.gz"


#path to the output directory for consensus sequence alignments. Make prior to running.
outdir="/work/bs66/dasanthera_novaseq/consensus_alignments"


#set up list of samples for which to generate consensus sequences
#this effectively lists sample names as listed in the .vcf, if following same protocol
#if sample names were changed at some point in .vcf, this will throw an error
samplelist=($mapped_filtered_reads/*.bam)
insample="${samplelist[$SLURM_ARRAY_TASK_ID]}"
########################


#1. filter gvcf
#2. index the reference for annotations of interest
#3. pipe this into bcftools consensus to generate sequence alignment for samples


###vcftools code -- sort that after I make sure the second part is working
#could also consider just using the GT plugin..? Might make more sense,

#indexing for annotations and generating consensus sequences
samtools faidx --region-file $annotations_file $refgenome | 
 bcftools consensus $invcf \
 --iupac-codes \
 --missing "N" \
 --mark-del "-" \
 -s $insample\
 -o $outdir/test_consensus.fa



