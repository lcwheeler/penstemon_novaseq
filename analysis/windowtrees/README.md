### Sliding window tree analysis

#### convert vcf files to fasta alignments
* install mvftools (probably better software to do this but I was trying something else)
* See `array_mvftools_vcf2fasta.sh`

#### split fasta alignments into 50kb stretches
* Perform this on scaffold-by-scaffold basis
* Enables missing data threshold: windows in which an individual has >= proportion of missing data (user-specified) 
* See `create_fasta_window_alignments.py`
    * Example usage:
    * `python3 create_fasta_window_alignments.py alignment.scaffold_1085.fa 100000 scaffold_1085`

#### estimate trees (and models of nucleotide substitution) for each 50kb stretch with iqtree
* Perform with array script
* P. montanus as outgroup
* See `array_iqtree_50kb_windows.sh`


#### make species tree in ASTRAL with 50kb window trees
* first, cat gene trees for each window into a single tree file
```shell
cd 50kb_alignments
for i in *;
do
    cat $i/*.treefile >> combined_50kb_trees.tre
done
```
* then, set up namefile, to assign taxa in gene trees to species tree
    * see `namemap_astral_v1.txt`
* run ASTRAL
    * see `run_astral_50kb_windows.sh`
    `java -jar $astral -i $treefile -a $outpath/namemap_astral_v1.txt --outgroup P_montanus -o $outpath/species_tree_astral.tre`
