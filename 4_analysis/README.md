## Analysis

### pixy
`cat *dxy.txt > allscaf_dxy.txt`
`cat *pi.txt > allscaf_pi.txt`
`cat *fst.txt > allscaf_fst.txt`


### Trees

#### Prepare fasta files for gene tree inference
The fasta files generated from consensus have newlines every 50bp. Additionally, they are sorted by species. We want to have a fasta file for each scaffold, with the sequence for that scaffold from each species. There is a script, [`TREES_1.rearrange_consensus_sequences.sh`](TREES_1.rearrange_consensus_sequences.sh), which does this. Note that this functions on Linux OS, and will need modified slightly if using on MacOS (read script for more details). Script will need edited to match variable naming scheme. Use:

`bash TREES_1.rearrange_consensus_sequences.sh`



#### ASTRAL species tree
##### Generate gene tree input files
This is done in two parts, with three files:
1. Make output directories for each scaffold. 
2. Generate gene trees with [`TREES_2.ARRAY_generate_genetree_infiles_masterscript.sh`](TREES_2.ARRAY_generate_genetree_infiles_masterscript.sh), which uses the custom python script [`TREES.create_fasta_window_alignments.py`](TREES.create_fasta_window_alignments.py). The missing data threshold is set in the python script; otherwise, it should function without modification. Parameters are defined in the batch script, including:
	* Input fasta file
	* Window size
	* prefix to append to outfiles (can function as outdir)

##### Estimate gene trees
This pipeline uses [IQtree](http://www.iqtree.org/) for gene tree inference. Again, this is done in two parts: 
1. Make output directories for gene trees estimated along each scaffold.
2. Estimate gene trees (in array batch submission), setting outgroup (here it is P. montanus), and specifying substitution model inference and out prefix.
* See [`TREES_3.ARRAY_estimate_genetrees_iqtree.sh`](TREES_3.ARRAY_estimate_genetrees_iqtree.sh)


##### Estimate species tree in ASTRAL
* cat trees into one large file. Change to directory with subdirectories, each containing gene trees estimated for each scaffold.

```shell
cd /work/bs66/dasanthera_novaseq/analysis/genetree_outfiles
for i in scaf_*;
do
    cat $i/*.treefile >> combined_windowtrees.tre
done
```

* also, keep a log of which trees were catted to the treefile, and in what order. This will come in handy later. It is a good idea to verify that the order is the same, but if done this way it should be correct.
```shell
for i in scaf_*;
do
    echo $i/*.treefile | tr " " "\n" >> tmpout.txt
done

cat --number tmpout.txt > numbered_treepaths.txt
rm tmpout.txt
```
* See Astral species tree script





