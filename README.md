## NovaSeq data for 18 accessions of *Penstemon* subgenus *Dasanthera*

### Note: this pipeline requires the use of conda installs of various softwares. For info on navigating installation of these packages with conda on the cluster at South Carolina, see [this readme file](https://github.com/benstemon/dasanthera_novaseq/blob/main/conda_info/README.md)

### Quality control pipeline
1. Merge Illumina lanes by forward and reverse reads
2. Run QC of raw sequencing reads (fastqc)
3. Trim adapters, quality filter, and enable base correction in overlapped regions (fastp)
4. Run QC on trimmed, filtered data (fastqc)
5. Summarize results (multiqc)

#### 1. Merge data across lanes
* See [`merge_illumina_lanes.sh`](https://github.com/benstemon/dasanthera_novaseq/blob/main/QC/merge_illumina_lanes.sh)

```shell
#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=merge_lanes

#Note: I found this script from here:
#https://www.biostars.org/p/317385/

cd $SLURM_SUBMIT_DIR


for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)

    do echo "Merging R1"

cat "$i"_L00*_R1_001.fastq.gz > "$i"_merged_L001_R1_001.fastq.gz

       echo "Merging R2"

cat "$i"_L00*_R2_001.fastq.gz > "$i"_merged_L001_R2_001.fastq.gz

done;
```

#### 2. QC on raw reads with fastqc

#### 3. Quality filtering with fastp
* Default options for quality filtering
* Base correction for overlapping reads enabled
* poly-x trimming on 3' ends enabled
* Unmerged reads that pass filter are included in separate files (default is to delete. reads failing to pass filter are not included)
* See [`run_fastp_novaseq.sh`](https://github.com/benstemon/dasanthera_novaseq/blob/main/QC/run_fastp_novaseq.sh)

```shell
#!/bin/sh

#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p wessinger-48core
#SBATCH --job-name=testrun_fastp

cd $SLURM_SUBMIT_DIR


#path to fastp
fastpdir="/work/bs66/software"

#path to illumina reads
reads="/work/bs66/dasanthera_novaseq/merged_reads"


#fastp for loop
for r1in in $reads/*_R1_001.fastq.gz; 
do
    r2in="${r1in/R1_001.fastq.gz/R2_001.fastq.gz}"
    r1out="${r1in##*/}"
    r2out="${r1out/R1_001.fastq.gz/R2_001.fastq.gz}"
    $fastpdir/./fastp -i "$r1in" -I "$r2in" --out1 "${r1out/merged_L001_R1_001.fastq.gz/trimmed_L001_R1_001.fastq.gz}" --out2 "${r1out/merged_L001_R1_001.fastq.gz/trimmed_L001_R2_001.fastq.gz}" --unpaired1 "${r1out/merged_L001_R1_001.fastq.gz/unpaired_L001_R1_001.fastq.gz}" --unpaired2 "${r1out/merged_L001_R1_001.fastq.gz/unpaired_L001_R2_001.fastq.gz}" -x -c -w 16 -h "${r1out/merged_L001_R1_001.fastq.gz/html}" -j "${r1out/merged_L001_R1_001.fastq.gz/json}"
done
```

This produces trimmed reads, that have been paired (in *trimmed*.fastq.gz) and unpaired (in *unpaired*.fastq.gz).
Moving forward with mapping, these unpaired reads can be mapped in the same way as paired reads, and then downstream BAMs can be merged with Picard's MergeSamFiles (or some other approach)

#### 4. fastqc on trimmed reads
* see [`run_fastqc_novaseq_trimmed.sh`](https://github.com/benstemon/dasanthera_novaseq/blob/main/QC/run_fastqc_novaseq_trimmed.sh)

```shell
#!/bin/sh

#SBATCH -N 1
#SBATCH -n 20
#SBATCH -p wessinger-48core
#SBATCH --job-name=fastqc_trimmed

cd $SLURM_SUBMIT_DIR

#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile

#activate conda environment with fastqc and multiqc installed
conda activate QC


#path to filtered illumina reads
reads="/work/bs66/dasanthera_novaseq/filtered_reads"
outdir="/work/bs66/dasanthera_novaseq/fastqc_filtered"

#move to directory with reads
cd $reads

#store fastq files as array
files=(*.fastq.gz)

#perform fastqc -- distinction from for loop is this can process -t files simultaneously
fastqc "${files[@]}" -t 20 -o $outdir
```
