# Step by step tutorial based on the h3achipimputation pipeline

The actual imputation tutorial describes steps to follow to perform a genotype imputation using minimac4 based on the h3achipimputation pipeline on https://github.com/h3abionet/chipimputation.  
All consecutive steps (commands given in 'cmd COMMAND' sections) must be run to ensure high-quality results.

## 1. Requirements

Step 1 defines the requirements for this tutorial (e.g. required software packages and reference genome files) and suggests example commands on how to process the files into suitable formats.

### 1.1.  Software packages

Required software packages for this tutorial are listed below with the versions used. 

BCFtools v1.9 [https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2](https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2)  
R v3.4.1 [https://www.r-project.org/](https://www.r-project.org/).  
R package data.table [https://github.com/Rdatatable/data.table/wiki/Installation](https://github.com/Rdatatable/data.table/wiki/Installation)  
R package sm [https://cran.r-project.org/package=sm](https://cran.r-project.org/package=sm)  
Eagle v2.4 [https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.tar.gz](https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.tar.gz)  
Minimac4 [https://github.com/statgen/Minimac4](https://github.com/statgen/Minimac4)

However, a singularity container with all required tools will be used.
The singularity image can be downloaded from [shub://h3abionet/chipimputation:minimac4](shub://h3abionet/chipimputation:minimac4).  
Ensure you have singularity installed on your computer. This can be done as explained [here](https://www.sylabs.io/docs/).  

You will need to download and transfer the Singularity image:

```bash
singularity pull --name h3abionet-chipimputation-minimac4.simg shub://h3abionet/chipimputation
```

Once transferred, start your singularity container:

```bash
singularity shell -B /data:data  h3abionet-chipimputation-minimac4.simg
```

// TODO add command explanations

### 1.2. Reference genome and genetic map files
 
#### 1.2.1 Fasta files
Homo Sapiens assembly hg37 is used and the required files are `human_g1k_v37_decoy.fasta` and `human_g1k_v37_decoy.fasta.fai`  
The files are available for [downloading](https://software.broadinstitute.org/gatk/download/bundle) at Broad Insitute FTP server at [ftp://ftp.broadinstitute.org/bundle/b37/](ftp://ftp.broadinstitute.org/bundle/b37/).
 
#### 1.2.2. Genetic map files for phasing with Eagle
Genetic map file (all chromosomes in a single file) with recombination frequencies is downloaded together with Eagle for GRCh37/hg19.  
`genetic_map_hg19_withX.txt.gz`
The genetic map file is available for downloading at Eagle download page at 
[https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/](https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/).

### 1.3. Phasing and imputation reference panel files

#### 1.3.1 Obtain the reference panel files
For increased phasing and imputation accuracy, we recommend using a population-specific imputation reference panel, if available.  
If population-specific reference data is not available, for instance 1000 Genomes Project [(www.nature.com/articles/nature15393)](www.nature.com/articles/nature15393) data can be used instead.

For this tutorial, the 1000 Genomes Project data will be used as reference panel. 
`GRCh37/hg19` files are available at Minimac4 site already processed to be compatible with minimac4:
[https://genome.sph.umich.edu/wiki/Minimac4](https://genome.sph.umich.edu/wiki/Minimac4).
 - VCF files for phasing can be downloaded with the command 
    ```bash
    wget -np ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/*
    ``` 
- M3VCF files for imputation can be dowmloaded with the command 
    ```bash
    wget -np ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/*
    ``` 
  
 Save the phased reference panel VCF files per chromosome as follows (for consistency with the commands in the protocol):
 panel_phased_chr#.vcf.gz (where # is chromosome number)