#!/bin/sh

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p wessinger-48core
#SBATCH --job-name=variants_bcftools

###NOTE: this is set up as a batch array script. SLURM submit array code begins on line 50.


cd $SLURM_SUBMIT_DIR

#activate conda environment with packages installed
#needs samtools (v1.15.1 works) and bcftools (v1.15.1 works)
conda activate mapping_etc

#load vcftools (v0.1.17 on cluster works)
module load vcftools


#######
#SETUP#
#######

########################
#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile


#hard path to the reference genome and annotations file of interest
#mapped filtered reads should not have a "/" at the end. This should be the same directory as output from script: mapping_pipeline_array_s1a-f.sh
refgenome="/work/bs66/project_compare_genomes/annot_Pdavidsonii_genome.1mb.fasta"
annotations_file="/work/bs66/project_compare_genomes/annot_Pdavidsonii_gffread.genes.bed"
mapped_filtered_reads="/work/bs66/dasanthera_novaseq/mapped_filtered_bams"


#Location of filtered vcf, from script #3 in this pipeline.
invcf="/work/bs66/dasanthera_novaseq/GVCF_VCFs/CHANGETHIS.vcf.gz"


#path to the output directory for consensus sequence alignments. Make prior to running.
outdir="/work/bs66/dasanthera_novaseq/consensus_alignments"


#number of cores used
numthreads=8


#set up list of samples for which to generate consensus sequences
#this effectively lists sample names as listed in the .vcf, if following same protocol
#if sample names were changed at some point in .vcf, this will throw an error
#infile is the full path to each input bam file
#insample is just the sample name (no path)
filelist=($mapped_filtered_reads/*.bam)
infile="${filelist[$SLURM_ARRAY_TASK_ID]}"

insample=${infile%_mapped_filtered.bam}
insample=${insample##*/}


########################







#2. index the reference for annotations of interest (not needed for full genome)
#3. pipe this into bcftools consensus to generate sequence alignment for samples




###vcftools code -- sort that after I make sure the second part is working


#generate consensus sequence
#if we want to pipe in 
bcftools consensus -f $refgenome \
 --iupac-codes \
 --absent "N" \
 --missing "N" \
 --mark-del "-" \
 --sample $infile \
 $invcf -o $outdir/test_consensus_fullgenome_onesample.fa



