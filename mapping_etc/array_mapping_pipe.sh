#!/bin/sh

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p wessinger-48core
#SBATCH --job-name=mapping_piped

#NOTE: this is an arrayed batch script. It requires special syntax to submit job.
#I have 18 total samples, and because counting starts at 0:
#submit with sbatch --array [0-17] array_mapping_pipe.sh


cd $SLURM_SUBMIT_DIR


#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

#activate conda environment with packages installed
#needs samtools v1.15.1 and bamutil v1.0.15 (I also installed bwa but that should be OK from module)
#version of samtools installed on cluster is old (no markdup and old fixmate options)

conda activate mapping_etc


#hard path to the filtered reads, reference genome, and out directory for mapped files
filtered_reads="/work/bs66/dasanthera_novaseq/filtered_reads"
refgenome="/work/bs66/project_compare_genomes/filtered_davidsonii_genome_1mb.fasta"
outdir="/work/bs66/dasanthera_novaseq/mapped_marked_bams"
statsdir="/work/bs66/dasanthera_novaseq/sumstats_marked_bams"


#identify target files
cd $filtered_reads

r1s=(*R1_001.fastq.gz)
read1="${r1s[$SLURM_ARRAY_TASK_ID]}"
read2="${read1/L001_R1/L001_R2}"

#pipeline:
#1. align
#2. convert sam to bam and fix read pairing info
#3. sort BAM by coordinates
#4. mark and remove duplicates
#5. index files

#mapping and cleaning pipeline
bwa mem -t 8 -M $refgenome $read1 $read2 | samtools fixmate -@ 8 -m -u -O bam - - | samtools sort -@ 8 -u | samtools markdup -@ 8 -u -s - - | bam clipOverlap --in -.ubam --out $outdir/"${read1/trimmed_L001_R1_001.fastq.gz/mapped_marked.bam}" --stats


#index reads
samtools index -b $outdir/"${read1/trimmed_L001_R1_001.fastq.gz/mapped_marked.bam}"



#COVERAGE
cd $outdir
#get information on coverage, with ASCII histogram
samtools coverage -m -A -w 40 "${read1/trimmed_L001_R1_001.fastq.gz/mapped_marked.bam}" > $statsdir/"${read1/trimmed_L001_R1_001.fastq.gz/coverage.txt}"

#you can tweak this to see how filtering out low quality mapped reads and duplicates would affect results
samtools coverage -m -A -w 40 -q 20 --ff "${read1/trimmed_L001_R1_001.fastq.gz/mapped_marked.bam}" > $statsdir/"${read1/trimmed_L001_R1_001.fastq.gz/coverage_MQ20_dedup.txt}"



#OTHER SUMMARY STATS
#piping just the summary numbers part of this command -- MQ 20
samtools stats -@ 8 -r $refgenome "${read1/trimmed_L001_R1_001.fastq.gz/mapped_marked.bam}" | grep -n "^SN" | cut -f 2- > $statsdir/"${read1/trimmed_L001_R1_001.fastq.gz/summarystats.txt}"


#filter for MQ 20 and no duplicates
samtools stats -@ 8 -d -q 20 -r $refgenome "${read1/trimmed_L001_R1_001.fastq.gz/mapped_marked.bam}" | grep -n "^SN" | cut -f 2- > $statsdir/"${read1/trimmed_L001_R1_001.fastq.gz/summarystats_MQ20_dedup.txt}"

