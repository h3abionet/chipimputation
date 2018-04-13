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
- The following executables should be available in PATH
  - [IMPUTE2](http://mathgen.stats.ox.ac.uk/impute/impute_v2.html) as `impute2`
  - [PLINK 1.9+](https://www.cog-genomics.org/plink2) as `plink2`
  - [VCFtools](https://vcftools.github.io/index.html) as `vcftools`
  - [BCFtools](https://samtools.github.io/bcftools/bcftools.html) as `bcftools`
  - [Eagle](https://data.broadinstitute.org/alkesgroup/Eagle/) as `eagle`

## Getting started
 1. Clone the repo
 2. Run the test imputation
```
nextflow run imputation.nf -C imputation_nf.test.config
```
 3. check for results in `outfolder`
