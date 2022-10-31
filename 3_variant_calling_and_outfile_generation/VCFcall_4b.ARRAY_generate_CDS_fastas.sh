#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p defq
#SBATCH --job-name=pullCDS
#SBATCH --output=slurm-generate-CDS-fasta_%j.out


#This is an arrayed batch script, requiring special syntax to submit. The info needed to determine array number begins on line 27. Submit with sbatch --array [0-n-1], where n is the number of consensus genome fastas in your input directory (genomedir). For me this is sbatch --array [0-17]


cd $SLURM_SUBMIT_DIR

#load bedtools
module load bedtools


#path to CDS.bed file from the annotation, and desired outfile
bedfile="/work/bs66/project_compare_genomes/annot_Pdavidsonii_gffread.genes.bed"
outdir="/work/bs66/dasanthera_novaseq/consensus_alignments/CDS_fastas"


#path to directory with genomes (fasta format) for each sample
genomedir="/work/bs66/dasanthera_novaseq/consensus_alignments/individual_fullgenome_fastas"


#set up variables array
filelist=($genomedir/*.fa)
infile="${filelist[$SLURM_ARRAY_TASK_ID]}"
samplename="${infile##*/}"
samplename="${samplename/consensus_fullgenome_/}"



#bedtools command to extract CDS for all mapped individuals
bedtools getfasta -fi $infile -bed $bedfile -fo $outdir/CDS_$samplename

