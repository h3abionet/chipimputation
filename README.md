# Imputation Workflow

This repo contains the workflow which was developed as part of
the H3ABioNet Hackathon held in Pretoria, SA in 2016.
 - We track our open tasks using github's [issues](https://github.com/h3abionet/chipimputation/issues)
 - The 1000ft view is located on [our Trello board](https://trello.com/b/Dp08chq7/stream-d-imputation-and-phasing).

## Setup (native cluster)

#### Headnode
  - [Nextflow](https://www.nextflow.io/) (can be installed as local user)
   - NXF_HOME needs to be set, and must be in the PATH
   - Note that we've experienced problems running Nextflow when NXF_HOME is on an NFS mount.
   - The Nextflow script also needs to be invoked in a non-NFS folder
  - Java 1.8+

#### Compute nodes

- The compute nodes need access to shared storage for input, references, output
- The following commands need to be available in PATH on the compute nodes

  - `impute2` from [IMPUTE2](http://mathgen.stats.ox.ac.uk/impute/impute_v2.html)
  - `plink` from [PLINK 1.9+](https://www.cog-genomics.org/plink2)
  - `vcftools` from [VCFtools](https://vcftools.github.io/index.html)
  - `bcftools`from [bcftools](https://samtools.github.io/bcftools/bcftools.html)
  - `bgzip` from [htslib](http://www.htslib.org)
  - `eagle` from [Eagle](https://data.broadinstitute.org/alkesgroup/Eagle/)
  - `python2.7`

## Getting started

### Basic test (testing reading / writing to shared home dir)
 1. Clone this repo and switch to this branch
```
cd $HOME 
https://github.com/h3abionet/chipimputation.git 
cd chipimputation
git branch azure-test
```

 2. Run the "tiny" dataset included
```
nextflow -log $HOME/chipimputation/nextflow.log run $HOME/chipimputation/imputation.nf -c $HOME/chipimputation/nextflow.test.tiny.config -w $HOME/chipimputation/nf-workdir
```
 3. check for results in `outfolder`
```
wc -l $HOME/chipimputation/output/impute_results/FINAL_VCFS/*
```

### Larger dataset (reading / writing to shared datastore that allows for fast reading and writing) 
 1. Download this slightly larger dataset: [small.tar.bz2](https://goo.gl/cYk51U) and extract into the `samples` folder
 2. Edit `dir` and `studyDir` in `nextflow.test.tiny.config` to point to the correct samples directory. 
 2. Run this "small" dataset with (change the file paths to your setup)
```
nextflow -log /global5/scratch/gerrit/projects/chipimputation/test.small/nextflow.log run $HOME/chipimputation/imputation.nf -c $HOME/chipimputation/nextflow.test.small.config -w /global5/scratch/gerrit/projects/chipimputation/test.small/nf-workdir -profile pbs -resume
```
 3. check for results in `outfolder`
```
wc -l /global5/scratch/gerrit/projects/chipimputation/test.small/nf-output/impute_results/FINAL_VCFS/*
```
