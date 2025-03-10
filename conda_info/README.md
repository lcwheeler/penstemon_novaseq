# Installing QC packages in conda

## Create QC env with python 3.7
```
module load python3/anaconda/2019.07
conda create --name QC
```

## Activate fresh env
`conda activate QC`

## Add appropriate channels to conda config
```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

## Install software, with newest version included.
I was having issues previously with older versions and/or incompatibilities with libgc version 2.17. Adding channels (I think especially conda-forge) resolved this issue for me.

```
conda install -y -c bioconda multiqc multiqc=1.6
conda install -y -c bioconda fastqc fastqc=0.11.8
conda install -y -c bioconda fastp fastp=0.23.2
```
