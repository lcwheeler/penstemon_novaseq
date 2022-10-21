#!/bin/sh

#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p wessinger-48core
#SBATCH --job-name=pixy_10kb

##NOTE: this is an arrayed batch script, requiring special syntax to submit jobs. The array in this case is relevant for the number of scaffolds (chromosomes) on which you would like to run pixy. Changes to the parameter "chromlist" are also necessary.


cd $SLURM_SUBMIT_DIR


#source bash profile and activate conda environment with pixy installed
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc


#paths to (1) the input vcf, (2) the desired out-directory, and (3) the populations file
invcf="/work/bs66/dasanthera_novaseq/GVCF_VCFs/filtered_consensus_ready.vcf.gz"
outdir="/work/bs66/dasanthera_novaseq/analysis/pixyout_fullgenome_10kb"
popfile="/work/bs66/dasanthera_novaseq/analysis/popfile_pixy.txt"


#List of scaffolds (chromosomes) to analyze
#An alternative (used here) is to read in the samtools faidx indexed file and pull scaffold names from there. Otherwise, list should be in the form ("chr1" "chr2" "chrN").
#When submitting the job, use sbatch --array [0-n], where n is the number of scaffolds -1
faidxfile="/work/bs66/project_compare_genomes/annot_Pdavidsonii_genome.1mb.fasta.fai"
chromlist=$(awk '{print $1}' $faidxfile)
chromlist=(${chromlist[@]})
#chromlist=("chr1" "chr2" "chr3" "chr4")
chromname="${chromlist[$SLURM_ARRAY_TASK_ID]}"


#Run pixy
pixy --stats pi fst dxy \
 --vcf $invcf \
 --populations $popfile \
 --chromosomes $chromname \
 --window_size 10000 \
 --n_cores 4 \
 --output_folder $outdir \
 --output_prefix scaffold_$chromname

