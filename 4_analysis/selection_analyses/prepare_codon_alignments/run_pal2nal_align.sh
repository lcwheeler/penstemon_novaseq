#!/bin/sh

#SBATCH -N 1
#SBATCH -n 12
#SBATCH -p wessinger-48core
#SBATCH --job-name=pal2nal
#SBATCH --output=slurm-pal2nal.out

# This script uses the pal2nal perl script to convert the CDS alignments to codon alignments
# Need to download the pal2nal script from http://www.bork.embl.de/pal2nal/ 

cd $SLURM_SUBMIT_DIR


#path to CDS.bed file from the annotation, and desired outfile
alndir="/work/lw74/habro/consensus_alignments_WithOutgroups_ReducedBED/CDS_fastas/rearranged_consensus_seqs_CDS_full"
pepdir="/work/lw74/habro/selection_analyses/CDS_reduced_BED_pep_aln"

mkdir -p "/work/lw74/habro/selection_analyses/CDS_reduced_BED_pal2nal_alignments"
outdir="/work/lw74/habro/selection_analyses/CDS_reduced_BED_pal2nal_alignments"

# path to the pal2nal perl script
pal2nal="/work/lw74/habro/selection_analyses/scripts/pal2nal.pl"


# iterate over the CDS fasta files and translate to AA seqs

#cd $pepdir
for f in $pepdir/*_AA.fasta
do
	name=$(basename $f)
	name=$(echo $name | sed 's/.rearranged.fa.relabel_AA.fasta//g')
	perl $pal2nal $f $alndir/$name".rearranged.fa.relabel.fasta" -output fasta -nogap -nomismatch > $outdir/$name".pal2nal"
done


