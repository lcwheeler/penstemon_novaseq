## NovaSeq data for 18 accessions of *Penstemon* subgenus *Dasanthera*

### Note: this pipeline requires the use of conda installs of various softwares. For info on navigating installation of these packages with conda on the cluster at South Carolina, see [this readme file](https://github.com/benstemon/dasanthera_novaseq/blob/main/conda_info/README.md)

### Quality control pipeline
1. Run QC of raw sequencing reads (fastqc)
2. Merge Illumina lanes by forward and reverse reads
4. Trim adapters, quality filter, and enable base correction in overlapped regions (fastp)
5. Run QC on trimmed, filtered data (fastqc)
6. Summarize results (multiqc)



#### 1. QC on raw reads with fastqc, summarized with multiqc
* see [`run_QC_s1.sh`](https://github.com/benstemon/dasanthera_novaseq/blob/main/QC/run_QC_s1.sh)


#### 2. Merge data across lanes
* See [`merge_illumina_lanes.sh`](https://github.com/benstemon/dasanthera_novaseq/blob/main/QC/merge_illumina_lanes.sh)





#### 3-5. Quality filtering with fastp, and check QC of trimmed reads
* Default options for quality filtering
* Base correction for overlapping reads enabled
* poly-x trimming on 3' ends enabled
* Unmerged reads that pass filter are included in separate files (default is to delete. reads failing to pass filter are not included)
* See [`run_QC_s3-5.sh`](https://github.com/benstemon/dasanthera_novaseq/blob/main/QC/run_QC_s3-5.sh)


This produces trimmed reads, that have been paired (in *trimmed*.fastq.gz) and unpaired (in *unpaired*.fastq.gz).
Moving forward with mapping, these unpaired reads can be mapped in the same way as paired reads, and then downstream BAMs can be merged with Picard's MergeSamFiles (or some other approach) 
