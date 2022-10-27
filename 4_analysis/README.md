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
2. Generate gene trees with [`TREES_2.ARRAY_generate_genetree_infiles_masterscript.sh`](TREES_2.ARRAY_generate_genetree_infiles_masterscript.sh), which uses the custom python script [`TREES.create_fasta_window_alignments.py`](TREES.create_fasta_window_alignments.py). This should function without modification. Parameters are defined in the batch script, including:
	* Input fasta file
	* Window size
	* Missing data threshold (only generates windows in which all species pass missing data threshold)
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
    cat $i/*.treefile >> combined_10kbwindowtrees.tre
done
```

* also, keep a log of which trees were catted to the treefile, and in what order. This will come in handy later. It is a good idea to verify that the order is the same, but if done this way it should be correct.
```shell
for i in scaf_*;
do
    echo $i/*.treefile | tr " " "\n" >> tmpout.txt
done

cat --number tmpout.txt > numbered_10kbtreepaths.txt
rm tmpout.txt
```

I put both of these outputs in a new directory, treemetrics. They will be used for several downstream analyses.


* Estimate the species tree in Astral. See [`SPECIESTREE_1.astral_windowtrees.sh`](SPECIESTREE_1.astral_windowtrees.sh)


#### Tree metrics
##### Robinson-Foulds distance
Using the conda install version of ete3 here. I am having trouble getting this to work through batch submission, so used an interactive node instead.

```
idev

#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile


#activate conda environment with packages installed
#need ete3 installed
conda activate quibl_etc


#specify reference (species) tree, and the combined window tree file
#reference tree should not have annotations so I made a new, cleaned up treefile
reftree="/work/bs66/dasanthera_novaseq/analysis/astral_trees/astral_10kb_noannotations.tre"
treelist="/work/bs66/dasanthera_novaseq/analysis/treemetrics/combined_10kbwindowtrees.tre"
outdir="/work/bs66/dasanthera_novaseq/analysis/treemetrics"


#run ete3
ete3 compare --src_tree_list $treelist -r $reftree --unrooted --taboutput > $outdir/RFdistance_10kbwindow_astralref.txt
```


