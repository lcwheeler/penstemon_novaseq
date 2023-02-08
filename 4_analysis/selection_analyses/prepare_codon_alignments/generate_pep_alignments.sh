#!/bin/sh

#SBATCH -N 1
#SBATCH -n 12
#SBATCH -p wessinger-48core
#SBATCH --job-name=aa_convert
#SBATCH --output=slurm-generate-aa-convert.out

# This script uses the MACSE aligner to convert the gapless CDS alignments extracted from full genome VCF to AA/pep sequence alignments


cd $SLURM_SUBMIT_DIR


#path to CDS.bed file from the annotation, and desired outfile
alndir="/work/lw74/habro/consensus_alignments_WithOutgroups_ReducedBED/CDS_fastas/rearranged_consensus_seqs_CDS_full"
mkdir -p "/work/lw74/habro/selection_analyses/CDS_reduced_BED_pep_aln"
outdir="/work/lw74/habro/selection_analyses/CDS_reduced_BED_pep_aln"


# iterate over the CDS fasta files and translate to AA seqs

cd $alndir

for f in *.relabel.fasta
do
	java -jar /work/lw74/software_installs/macse/macse_v2.06.jar -prog translateNT2AA -seq $f
done


mv *_AA.fasta $outdir

