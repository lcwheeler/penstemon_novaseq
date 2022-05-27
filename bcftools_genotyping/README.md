# Directory for genotype and SNP calling

Pipeline:
* Genotype likelihoods with bcftools mpileup
    * Contig-by-contig basis
    * Max depth of 100x
* Variants called with bcftools call
* Reads not passing quality filters (marked previously or ambiguous genotypes) filtered to missing data with bcftools filter
