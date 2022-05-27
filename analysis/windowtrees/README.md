### Sliding window tree analysis

#### convert vcf files to fasta alignments
* install mvftools (probably better software to do this but I was trying something else)
* See array_mvftools_vcf2fasta.sh

#### split fasta alignments into 100kb stretches
* Perform this on scaffold-by-scaffold basis
* Enables missing data threshold: windows in which an individual has >= proportion of missing data (user-specified) 
* See create_fasta_window_alignments.py

