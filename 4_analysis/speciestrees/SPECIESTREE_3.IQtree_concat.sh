#!/bin/sh

#SBATCH -N 1
#SBATCH -n 32
#SBATCH -p defq
#SBATCH --job-name=concat_tree
#SBATCH --output=iqtree_speciestree_%j.out


cd $SLURM_SUBMIT_DIR

#load conda env. Need numpy, matplotlib, biopython, (others?). installed.
source /home/lw74/.bashrc
source /home/lw74/.bash_profile


#path to scaffold fasta alignments
scaffolds="/work/lw74/habro/consensus_alignments_WithOutgroups/CDS_fastas/filtered_consensus_fasta/passed_files"


#path to output file
outpath="/work/lw74/habro/species_trees_habro"


#infer species tree in IQtree
iqtree -p $scaffolds \
 -m GTR+I+R \
 -o THD_hirs_005J_S373 \
 -T 32 \
 -B 1000 \
 --prefix $outpath/concat_speciestree_CDS

