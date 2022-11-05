## Analysis

### pixy
`cat *dxy.txt > allscaf_dxy.txt`
`cat *pi.txt > allscaf_pi.txt`
`cat *fst.txt > allscaf_fst.txt`


### Trees

#### Windowed gene trees
##### Prepare fasta files for windowed gene tree inference
The fasta files generated from consensus have newlines every 50bp. Additionally, they are sorted by species. We want to have a fasta file for each scaffold, with the sequence for that scaffold from each species. There is a script, [`TREES_1.rearrange_consensus_sequences.sh`](genetrees/TREES_1.rearrange_consensus_sequences.sh), which does this. Note that this functions on Linux OS, and will need modified slightly if using on MacOS (read script for more details). Script will need edited to match variable naming scheme. Use:

`bash TREES_1.rearrange_consensus_sequences.sh`


##### Generate windowed gene tree input files
This is done in two parts, with three files:
1. Make output directories for each scaffold. 
2. Generate gene trees with [`TREES_2.ARRAY_generate_windowed_genetree_infiles_masterscript.sh`](genetrees/TREES_2.ARRAY_generate_windowed_genetree_infiles_masterscript.sh), which uses the custom python script [`TREES.create_fasta_window_alignments.py`](TREES.create_fasta_window_alignments.py). This should function without modification. Parameters are defined in the batch script, including:
	* Input fasta file
	* Window size
	* Missing data threshold (only generates windows in which all species pass missing data threshold)
	* prefix to append to outfiles (can function as outdir)


##### Estimate windowed gene trees
This pipeline uses [IQtree](http://www.iqtree.org/) for gene tree inference. Again, this is done in two parts: 
1. Make output directories for gene trees estimated along each scaffold.
2. Estimate gene trees (in array batch submission), setting outgroup (here it is P. montanus), and specifying substitution model inference and out prefix.
* See [`TREES_3.ARRAY_estimate_windowed_genetrees_iqtree.sh`](genetrees/TREES_3.ARRAY_estimate_genetrees_iqtree.sh)



#### CDS gene trees

##### Prepare CDS fasta files for gene tree inference
The fasta files generated for CDS regions are initially sorted by sample. To estimate CDS gene trees, we need a fasta for each CDS, with each sample's sequence included. First we will want to reformat our fasta files, so they are single line fasta rather than interleaved. Use this code on the CDS fastas.
```shell

for i in *.fa;
do
    awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < $i > "${i/.fa/.fixed.fa}"
    rm $i
done

```

Next, we will want to generate the infiles for gene tree inference. To generate these, see [`CDS_TREES.concat_CDS_fastas.py`](genetrees/CDS_TREES.concat_CDS_fastas.py). This python script takes input fasta, output directory, and missing data threshold, and appends to an output fasta for each CDS, naming the output after the scaffold and region the CDS corresponds to. A simple shell for loop can be run with the python script to add all samples to the output. The missing data threshold here filters for individuals, rather than windows (i.e., if a sample has more missing data than desired, that individual is not added to the output fasta, rather than the output fasta not being generated).

```shell
for i in individual_CDS_fastas/CDS_*.fa;
do
    python3 CDS_TREES.concat_CDS_fastas.py -i $i -o /work/bs66/dasanthera_novaseq/analysis/CDS_genetree_infiles -m 0.5
done
```

Here I generate aligned fasta for each CDS, excluding individuals with > 50% missing data. The script is fast, and possible to use without batch submission. But for many samples or many CDS regions, it would be wise to write a batch script so as not to overload the head node.


##### Estimate CDS gene trees
Still using IQtree for gene tree inference here. This used to be done on a scaffold-by-scaffold basis, but because that information has been lost in the header, instead we are just submitting all files for inference. They could be split up into sections to speed this part up.
* To infer gene trees, see [`CDS_TREES_1.estimate_CDS_genetrees_iqtree.sh`](genetrees/CDS_TREES_1.estimate_CDS_genetrees_iqtree.sh)




#### ASTRAL species tree: windowed analysis

First, cat trees into one large file. Change to directory with subdirectories, each containing gene trees estimated for each scaffold.

```shell
cd /work/bs66/dasanthera_novaseq/analysis/genetree_outfiles
for i in scaf_*;
do
    cat $i/*.treefile >> combined_10kbwindowtrees.tre
done
```

Also, keep a log of which trees were catted to the treefile, and in what order. This will come in handy later. It is a good idea to verify that the order is the same, but if done this way it should be correct.

```shell
for i in scaf_*;
do
    echo $i/*.treefile | tr " " "\n" >> tmpout.txt
done

cat --number tmpout.txt > numbered_10kbtreepaths.txt
rm tmpout.txt
```

I put both of these outputs in a new directory, `treemetrics`. They will be used for several downstream analyses. Next, use the combined windowtrees file to estimate the species tree in ASTRAL.

* See [`SPECIESTREE_1.astral_windowtrees.sh`](speciestrees/SPECIESTREE_1.astral_windowtrees.sh)



#### ASTRAL species tree: CDS analysis

Perform the same commands as described for the windowed analysis; cat trees to a combined tree file and keep a log of the order in which this was done. 

```shell
cd /work/bs66/dasanthera_novaseq/analysis/CDS_genetree_outfiles
for i in *.treefile;
do
    cat $i >> combined_CDStrees.tre
done

for i in *.treefile;
do
    echo $i | tr " " "\n" >> tmpout.txt
done

cat --number tmpout.txt > numbered_CDStreepaths.txt
rm tmpout.txt

mv combined_CDStrees.tre ../treemetrics
mv numbered_CDStreepaths.txt ../treemetrics
```

Then, estimate the species tree in ASTRAL.
* See [`SPECIESTREE_2.astral_CDS.sh`](speciestrees/SPECIESTREE_2.astral_CDS.sh)




#### Concatenated ML species tree
The concatenated ML species tree is estimated in IQtree. I specified the GTR+I+R model, and estimated rates for each scaffold, performing 1000 ultra-fast bootstrap replicates, and specifying _P. montanus_ as the outgroup.
* See [`SPECIESTREE_3.IQtree_concat.sh`](speciestrees/SPECIESTREE_3.IQtree_concat.sh)


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


