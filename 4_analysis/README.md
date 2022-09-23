## Analysis

### pixy
`cat *dxy.txt > allscaf_dxy.txt`
`cat *pi.txt > allscaf_pi.txt`
`cat *fst.txt > allscaf_fst.txt`


### trees

#### Prepare fasta files for gene tree inference
The fasta files generated from consensus have newlines every 50bp. Additionally, they are sorted by species. We want to have a fasta file for each scaffold, with the sequence for that scaffold from each species. There is a script, [`rearrange_consensus_sequences.sh`](rearrange_consensus_sequences.sh), which does this.Note that this functions on Linux OS, and will need modified slightly if using on MacOS (read script for more details). Script will need edited to match variable naming scheme. Use:

`bash rearrange_consensus_sequences.sh`


#### Generate gene tree input files
* python script
* master bash script

