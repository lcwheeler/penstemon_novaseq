#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=convert_and_Dstat
#SBATCH --output=slurm-Dsuite_Dtrios_%j.out
#SBATCH --error=slurm-Dsuite_Dtrios_%j.error


###NOTE: this script can only be run once the first step (DSTATS_1) has been completed, as it relies on output from those analyses


cd $SLURM_SUBMIT_DIR


#######
#SETUP#
#######

#source appropriate environments to enable use of conda installs through batch submission
source /home/lw74/.bashrc
source /home/lw74/.bash_profile


#set Dsuite variable to path, and for f-branch test, the path to directory with dtools.py
dtools="/work/lw74/software_installs/Dsuite/utils"


#specify invcfdir, invcf, outdir, faidx file, tree file
invcfdir="/work/lw74/habro/VCFs_WithOutgroups"
invcf="/work/lw74/habro/VCFs_WithOutgroups/filtered_consensus_ready_no-indels.vcf.gz"
dstat_outdir="/work/lw74/habro/Dstats_analyses/dsuite_outputs"
faidxfile="/work/lw74/refGenome/Pbar.2022.LG.fa.fai"
treefile="/work/lw74/habro/Dstats_analyses/CDS_WithOutgroups_intree_dtrios.pruned.relabel.tidy.tre"
test_trios="/work/lw74/habro/Dstats_analyses/test_trios.txt" 
popset="/work/lw74/habro/Dstats_analyses/popset_dtrios.txt"


##########
#ANALYSIS#
##########


#make array for specifying scaffolds -- used to specify input files, assuming they were generated on a per-scaffold basis
cd $dstat_outdir
chromlist=()
chromlist=$(awk '{print $1}' $faidxfile)
chromlist=(${chromlist[@]})
acount=${#chromlist[@]}


#for loop to append appropriate directories to the list
for((i=0;i<acount;i++));
do
    chromlist[i]="${chromlist[i]}".outfiles/"${chromlist[i]}"
done



#Combine results from Dtrios runs across scaffolds
Dsuite DtriosCombine -t $treefile -o WHOLEGENOME "${chromlist[@]}"



#do fbranch test in Dsuite for wholegenome
#note that this outputs only significant results at default p=0.01
#plots from the original paper put asterisks in significant blocks, so could
#change p-thresh to 1, and then put stars in blocks that are significant.
Dsuite Fbranch $treefile $dstat_outdir/WHOLEGENOME_combined_tree.txt > $dstat_outdir/WHOLEGENOME_Fbranch.txt




#move wholegenome files to new directories
mkdir WHOLEGENOME.outfiles
mv WHOLEGENOME_*.txt WHOLEGENOME.outfiles


#plot fbranch test output
cd WHOLEGENOME.outfiles
python $dtools/dtools.py WHOLEGENOME_Fbranch.txt $treefile
mv fbranch.png WHOLEGENOME.fbranch.png
mv fbranch.svg WHOLEGENOME.fbranch.svg
