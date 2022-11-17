#!/bin/sh

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p wessinger-48core
#SBATCH --job-name=mapping_pipeline

#NOTE: this is an arrayed batch script. It requires special syntax to submit job.
#I have 38 total samples, and because counting starts at 0:
#submit with sbatch --array [0-37] mapping_pipeline_array_s1a-f.sh


cd $SLURM_SUBMIT_DIR


#source appropriate environments to enable use of conda installs through batch submission
source /home/lw74/.bashrc
source /home/lw74/.bash_profile

#activate conda environment with packages installed
#needs samtools v1.15.1 and bamutil v1.0.15 (I also installed bwa but that should be OK from module)
#version of samtools installed on cluster is old (no markdup and old fixmate options)

#conda activate mapping_etc


#hard path to the filtered reads and reference genome
filtered_reads="/work/lw74/habro/filtered_reads"
refgenome="/work/lw74/refGenome/Pbar.2022.LG.fa"

#path to the out directories for filtered bams and summary stats. Make prior to running
mkdir -p /work/lw74/habro/mapped_filtered_bams
mkdir -p /work/lw74/habro/sumstats_mapped_filtered_bams

outdir="/work/lw74/habro/mapped_filtered_bams"
statsdir="/work/lw74/habro/sumstats_mapped_filtered_bams"
numthreads=8 #change this to match number of cores allocated


#identify target files
cd $filtered_reads

r1s=(*R1_001.fastq.gz)
read1="${r1s[$SLURM_ARRAY_TASK_ID]}"
read2="${read1/L001_R1/L001_R2}"


#mapping and cleaning pipeline:
#a. map with bwa mem
#b. fixmate -m fills in mate coords and add score tags to tell markdup which reads to keep
#c. sort alignment (put mapped reads in physical order)
#d. markdup marks and removes (-r) duplicate reads
#e. view filters to retain only reads with mapping quality >20 (99% mapping confidence)
#f. bam clipOverlap softclips read overlaps

bwa mem -t $numthreads -M $refgenome $read1 $read2 | \
 samtools fixmate -@ $numthreads -m -u -O bam - - | \
 samtools sort -@ $numthreads -u | \
 samtools markdup -r -@ $numthreads -u -s - - | \
 samtools view -@ $numthreads -h -q 20 -u | \
 bam clipOverlap --in -.ubam \
 --out $outdir/"${read1/trimmed_L001_R1_001.fastq.gz/mapped_filtered.bam}" --stats

#index reads
samtools index -b $outdir/"${read1/trimmed_L001_R1_001.fastq.gz/mapped_filtered.bam}"


#COVERAGE
cd $outdir
#get information on coverage, with ASCII histogram
samtools coverage -m -A -w 40 "${read1/trimmed_L001_R1_001.fastq.gz/mapped_filtered.bam}" > $statsdir/"${read1/trimmed_L001_R1_001.fastq.gz/coverage.txt}"

#OTHER SUMMARY STATS
#piping just the summary numbers part of this command 
samtools stats -@ $numthreads -r $refgenome "${read1/trimmed_L001_R1_001.fastq.gz/mapped_filtered.bam}" | \
 grep -n "^SN" | \
 cut -f 2- > $statsdir/"${read1/trimmed_L001_R1_001.fastq.gz/summarystats.txt}"

