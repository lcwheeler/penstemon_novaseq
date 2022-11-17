#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=pullCDS
#SBATCH --output=slurm-generate-CDS-fasta.out


#This is an arrayed batch script, requiring special syntax to submit. The info needed to determine array number begins on line 27. 
#Submit with sbatch --array [0-n-1], where n is the number of consensus genome fastas in your input directory (genomedir). For me this is sbatch --array [0-41]


cd $SLURM_SUBMIT_DIR

#load bedtools
module load bedtools


#path to CDS.bed file from the annotation, and desired outfile
bedfile="/work/lw74/refGenome/M4_annotation_putative_function_domain_added_blast_tomato.genemodels.noseq.bed"
mkdir -p /work/lw74/habro/consensus_alignments_WithOutgroups/CDS_fastas
outdir="/work/lw74/habro/consensus_alignments_WithOutgroups/CDS_fastas"


#path to directory with genomes (fasta format) for each sample
genomedir="/work/lw74/habro/consensus_alignments_WithOutgroups/individual_fullgenome_fastas"


#set up variables array
filelist=($genomedir/*.fa)
infile="${filelist[$SLURM_ARRAY_TASK_ID]}"
samplename="${infile##*/}"
samplename="${samplename/consensus_fullgenome_/}"



#gffread command to extract CDS for all mapped individuals
gffread -x $outdir/CDS_full_$samplename -g $infile $bedfile
