#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=mvftools_windowtree

##NOTE: this is an arrayed batch script. It requires special syntax to submit jobs.
#I have 12 total scaffolds, and because counting starts at 0:
#submit with sbatch --array [0-11] array_mapping_pipe.sh


cd $SLURM_SUBMIT_DIR

#load conda env. Need numpy, matplotlib, biopython, (others?). installed.
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc

#paths to important files and directories
mvftools="/work/bs66/software/mvftools" #path to mvftools directory
vcfdir="/work/bs66/dasanthera_novaseq/called_variants_bcftools"


#identify input vcfs
cd $vcfdir
readlist=(*.vcf)
vcfread="${readlist[$SLURM_ARRAY_TASK_ID]}"

#prepare .mvf output name
mvfread="${vcfread/.vcf/.mvf}"
fastaread="${vcfread/.vcf/.fa}"



#########################
#USING MVFTOOLS
# NOTE: MUST RUN CODE WITH PYTHON3: PYTHON 2 WILL NOT WORK. (my base python for this conda env is python3)
#########################

#convert the vcf to fasta format
python $mvftools/mvftools.py ConvertVCF2MVF --vcf $vcfread --out $mvfread
python $mvftools/mvftools.py ConvertMVF2Fasta --mvf $mvfread --out $fastaread

