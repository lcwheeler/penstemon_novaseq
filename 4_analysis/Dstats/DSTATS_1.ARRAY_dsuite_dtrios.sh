#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p defq
#SBATCH --job-name=convert_and_Dstat
#SBATCH --output=slurm-Dsuite_Dtrios_%j.out
#SBATCH --error=slurm-Dsuite_Dtrios_%j.error


###NOTE: this is set up as a batch array script. The array number is the number of scaffold.fasta files, -1
##0-11 for davidsonii 1mb


cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

#source appropriate environments to enable use of conda installs through batch submission
source /home/lw74/.bashrc
source /home/lw74/.bash_profile

#set Dsuite variable to path, and for f-branch test, the path to directory with dtools.py
dtools="/work/lw74/software_installs/Dsuite/utils"

#set Dsuite variable to path, and specify Dsuite popset file
popset="/work/lw74/habro/Dstats_analyses/popset_dtrios.txt"

#create a new output sub-directory
mkdir -p "/work/lw74/habro/Dstats_analyses/dsuite_outputs"

#specify invcfdir, invcf, outdir, faidx file, tree file
invcfdir="/work/lw74/habro/VCFs_WithOutgroups"
invcf="/work/lw74/habro/VCFs_WithOutgroups/filtered_consensus_ready_no-indels.vcf.gz"
dstat_outdir="/work/lw74/habro/Dstats_analyses/dsuite_outputs"
faidxfile="/work/lw74/refGenome/Pbar.2022.LG.fa.fai"
treefile="/work/lw74/habro/Dstats_analyses/CDS_WithOutgroups_intree_dtrios.pruned.relabel.tidy.tre"
test_trios="/work/lw74/habro/Dstats_analyses/test_trios.txt" 


##########
#ANALYSIS#
##########


#make array for specifying scaffolds -- to pull out appropriate scaffold, given array
chromlist=$(awk '{print $1}' $faidxfile)
chromlist=(${chromlist[@]})
chromname="${chromlist[$SLURM_ARRAY_TASK_ID]}"


#use bcftools to make scaffold-specific .vcf file
#bcftools view $invcf --regions $chromname -o $invcfdir/$chromname.vcf.gz -Oz

#index the vcf
#tabix $invcfdir/$chromname.vcf.gz


#Use Dsuite to estimate Dstatistics between all triplet pairs
#specifying the output prefix, the input vcf, and the popset
Dsuite Dtrios -o $dstat_outdir/$chromname -t $treefile $invcfdir/$chromname.vcf.gz $popset


#do fbranch test in Dsuite
#note that this outputs only significant results at default p=0.01
#plots from the original paper put asterisks in significant blocks, so could
#change p-thresh to 1, and then put stars in blocks that are significant.
Dsuite Fbranch $treefile $dstat_outdir/"${chromname}_tree.txt" > $dstat_outdir/"${chromname}_Fbranch.txt"

#move files to new directories
cd $dstat_outdir
mkdir $chromname.outfiles
mv "${chromname}"_* $chromname.outfiles

#plot fbranch test output
cd $chromname.outfiles
python $dtools/dtools.py "${chromname}_Fbranch.txt" $treefile
mv fbranch.png $chromname.fbranch.png
mv fbranch.svg $chromname.fbranch.svg

