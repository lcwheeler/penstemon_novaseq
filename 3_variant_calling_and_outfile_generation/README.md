## Genotyping and variant calling


### Prepare .bed file with annotated regions of interest
In this case, the .bed file produces genic regions, which we will use downstream to generate consensus sequences for genes of interest

* See [`VCFcall_1.gffread_convert_gff-to-bed.sh`](VCFcall_1.gffread_convert_gff-to-bed.sh)
* Requires installation of [gffread](https://github.com/gpertea/gffread)



### Call genotypes to produce a gVCF
The gVCF includes variants in addition to invariant sites (which are called in contiguous haplotype blocks). After this script has concluded (and the successful creation of the unfiltered VCF) you can remove the gVCF.

* see [`VCFcall_2.call_genotypes_mpileup_GVCF.sh`](VCFcall_2.call_genotypes_mpileup_GVCF.sh)
* Currently filters on base quality <20 and mapping quality <20, and calls invariant sites in such a way that you can filter on read depth and other parameters for invariant sites if you choose.



### Filter vcfs (goal: generate consensus sequence)
These filters are for the production of consensus sequences for use in a phylogenetic context. The filters used here are therefore appropriate for this purpose, but are likely not suitable for generating VCFs in a population genomic context.

Invariant and variant sites are filtered differently, because we want stricter filters on sites which are putatively variable. This means there are three main steps: (1) filtering invariant sites, (2) filtering variant sites, and (3) recombining post-filtered data into a single VCF. For this reason there are four provided scripts:

* [`VCFcall_3.filter_vcf_consensus.sh`](VCFcall_3.filter_vcf_consensus.sh) Performs all three steps in a single script. But it is slower, because it does not concurrently perform filtering steps.
* [`VCFcall_3a.filter_vcf_invariants_consensus.sh`](VCFcall_3a.filter_vcf_invariants_consensus.sh) Filters invariant sites. Currently no filters are actually placed on invariant sites, but one could reasonably try depth filters.
* [`VCFcall_3b.filter_vcf_variants_consensus.sh`](VCFcall_3b.filter_vcf_variants_consensus.sh) Filters variant sites. Current filters implemented include genotype quality <20, minimum mean depth 3x, maximum mean depth 40x, minimum depth (for individual genotype call) 2x, maximum missing data allowed = 20%. In addition, genotypes at sites with low quality variants (QUAL <20) are changed to match the reference allele.
* [`VCFcall_3c.merge_filtered_vcfs_consensus.sh`](VCFcall_3c.merge_filtered_vcfs_consensus.sh) Combines the two separate filtered VCFs into one large consensus VCF. From here one can move on to the production of consensus sequence alignments (fasta files) or perform certain analyses straight from the consensus VCF.



### Generate consensus sequence
This script uses the reference genome in conjunction with the consensus VCF file to generate consensus sequences for each sample. The command `samtools faidx $referencegenome` will need to be run prior to the use of this script. This is an ARRAY script, meaning that to run, you will need to implement some form of `sbatch --array [0-n] scriptname.sh`, where n is the number of samples -1, and `scriptname.sh` is a generic name for the actual script being run.
* See [`VCFcall_4a.ARRAY_generate_consensus_fullgenome.sh`](`VCFcall_4a.ARRAY_generate_consensus_fullgenome.sh`)



### Bonus: working with indels
Scripts to perform variant calling and consensus sequence generation while including indels can be found in [`filter_leave_indels`](filter_leave_indels/).

These are still a work in progress, because the consensus sequence method used here will not produce genome fasta files that are aligned to one another. On top of this, the indel detection method in bcftools doesn't appear to work super well; looking into alternative indel calling methods, but right now this is not a primary concern.



