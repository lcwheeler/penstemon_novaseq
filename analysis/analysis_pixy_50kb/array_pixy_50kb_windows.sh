#!/bin/sh

#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p wessinger-48core
#SBATCH --job-name=pixy_array

##NOTE: this is an arrayed batch script. It requires special syntax to submit jobs.
#I have 12 total scaffolds, and because counting starts at 0:
#submit with sbatch --array [0-11] array_pixy_50kb_windows.sh

cd $SLURM_SUBMIT_DIR

source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc

#paths to important files and directories
vcfdir="/work/bs66/dasanthera_novaseq/called_variants_bcftools"
outdir="/work/bs66/dasanthera_novaseq/analysis_pixy_50kb"
popfile="/work/bs66/dasanthera_novaseq/analysis_pixy_50kb/popfile_pixy.txt"


#identify input vcfs and contig name
cd $vcfdir
readlist=(*.vcf.gz)
vcfread="${readlist[$SLURM_ARRAY_TASK_ID]}"
contigname="${vcfread//[!0-9]/}"

#example pixy code:
pixy --stats pi fst dxy \
 --vcf $vcfread \
 --populations $popfile \
 --window_size 50000 \
 --n_cores 4 \
 --output_folder $outdir \
 --output_prefix scaffold_$contigname

