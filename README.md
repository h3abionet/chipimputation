# Imputation Workflow h3abionet/chipimputation

[![Build Status](https://travis-ci.org/h3abionet/chipimputation.svg?branch=master)](https://travis-ci.org/h3abionet/chipimputation)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker%20registry-Quay.io-red)](https://quay.io/h3abionet_org/imputation_tools)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8B%20%20%E2%97%8B-orange)](https://fair-software.eu)

## Introduction
Imputation is likely to be run in the context of a GWAS, studying population structure, and admixture studies. It is computationally expensive in comparison to other GWAS steps.
The basic steps of the pipeline is described in the diagram below:

![chipimputation pipeline workflow diagram](https://www.h3abionet.org/images/workflows/snp_imputation_workflow.png)

* The workflow is developed using [![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/) and imputation performed using [Minimac4](https://genome.sph.umich.edu/wiki/Minimac4). 
* It identifies regions to be imputed on the basis of an input file in VCF format, split the regions into small chunks, phase each chunk using the phasing tool [Eagle2](https://data.broadinstitute.org/alkesgroup/Eagle/) and produces output in VCF format that can subsequently be used in a [GWAS](https://github.com/h3abionet/h3agwas) workflow.
* It also produce basic plots and reports of the imputation process including the imputation performance report, the imputation accuracy, the allele frequency of the imputed vs of the reference panel and other metrics.    

**This pipeline comes with docker/singularity containers making installation trivial and results highly reproducible.**



## Getting started

### Running the pipeline with test dataset
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub.
You can run the pipeline using test data hosted in github with singularity without have to install or change any parameters.

```
nextflow run h3abionet/chipimputation/main.nf -profile test,singularity
```

- `test` profile will download the testdata from [here](https://github.com/h3abionet/chipimputation_test_data/tree/master/testdata_imputation)
- `singularity` profile will download the singularity image from [quay registry](https://quay.io/h3abionet_org/imputation_tools)

Check for results in `./output`

In case you are running the outdated version, run the code below prior to executing the pipeline.

```
nextflow pull h3abionet/chipimputation
```

### Start running your own analysis

Copy the `test.config` file from the `conf` folder by doing `cp <conf dir>/test.config .` and edit it to suit the path to where your files are stored.

Once you have edited the config file, run the command below.

```bash
nextflow run h3abionet/chipimputation -c "name of your config file" -profile singularity
```

- `singularity` profile will download the singularity image from [quay registry](https://quay.io/h3abionet_org/imputation_tools)

Check for results in `./output`


## Documentation
The h3achipimputation pipeline comes with detailed documentation about the pipeline.
This is found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. [Pipeline configuration](docs/configuration/config_files.md)  
    2.1. [Configuration files](docs/configs.md)  
    2.2. [Software requirements](docs/soft_requirements.md)  
    2.3. [Other clusters](docs/other_clusters.md)  
3. [Running the pipeline with test data](docs/usage.md)
4. [Running the pipeline with your own config](docs/usage.md)
5. [Running on local machine or cluster](docs/other_clusters.md)
6. [Running docker and singularity](docs/soft_requirements.md)


## Support
We track our open tasks using github's [issues](https://github.com/h3abionet/chipimputation/issues)


## Citation
This  workflow which was developed as part of the H3ABioNet Hackathon held in Pretoria, SA in 2016. Should want to reference it, please use:  
>Baichoo S, Souilmi Y, Panji S, Botha G, Meintjes A, Hazelhurst S, Bendou H, Beste E, Mpangase PT, Souiai O, Alghali M, Yi L, O'Connor BD, Crusoe M, Armstrong D, Aron S, Joubert F, Ahmed AE, Mbiyavanga M, Heusden PV, Magosi LE, Zermeno J, Mainzer LS, Fadlelmola FM, Jongeneel CV, Mulder N. Developing reproducible bioinformatics analysis workflows for heterogeneous computing environments to support African genomics. BMC Bioinformatics. 2018 Nov 29;19(1):457. doi: 10.1186/s12859-018-2446-1. PubMed PMID: 30486782; PubMed Central PMCID: [PMC6264621](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6264621/).
