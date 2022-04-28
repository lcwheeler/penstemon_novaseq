### NovaSeq data for 18 accessions of *Penstemon* subgenus *Dasanthera*

#### Step 1: Quality control pipeline
1. Run QC of raw sequencing reads (fastqc)
2. Trim adapters, quality filter, and enable base correction in overlapped regions (fastp)
3. Run QC on trimmed, filtered data (fastqc)
4. Summarize results (multiqc)
