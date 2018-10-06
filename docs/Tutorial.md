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

For increased phasing and imputation accuracy, we recommend using a population-specific imputation reference panel, if available.  
If population-specific reference data is not available, for instance 1000 Genomes Project [(www.nature.com/articles/nature15393)](www.nature.com/articles/nature15393) data can be used instead.

For this tutorial, the 1000 Genomes Project data will be used as reference panel. 
`GRCh37/hg19` files are available at Minimac4 site already processed to be compatible with minimac4:
[https://genome.sph.umich.edu/wiki/Minimac4](https://genome.sph.umich.edu/wiki/Minimac4).
 - VCF files for phasing can be downloaded with the command 
    ```bash
    wget -np ftp://share.sph.umich.edu/minimac3/G1K_P3_VCF_Files.tar.gz
    ``` 
- M3VCF files for imputation can be dowmloaded with the command 
    ```bash
    wget -np ftp://share.sph.umich.edu/minimac3/G1K_P3_M3VCF_FILES_WITH_ESTIMATES.tar.gz*
    ``` 
    
These phased reference panel VCF and M3VCF files are in chromosome.  
**NOTE**: The reference panel files should contain non-missing, phased genotypes!


## 2. Target/Chip data preparation

### 2.0. Get target data

For this tutorial, only a region on the genome will be used to reduce running time.
Extract chunk 6:428721-20428720 for `sample.vcf.gz`
```bash
bcftools index --tbi -f sample.vcf.gz
bcftools view --regions 6:428721-20428720 -m2 -M2 -v snps sample.vcf.gz -Oz -o sample_chr6_428721_20428720 .vcf.gz
```
> `bcftools index` parameters:   
--tbi generate index for VCF file  
`bcftools view` parameters:  
--regions region to extract of format `chromosome:startPostion-endPosition  
-m2 minimum allele count   
-M2 maximum allele count   
-v type of variants to include 
-Oz compressed output  

### 2.1. Check and remove invalid chromosomes

Only autosomal chromosomes are kept for this imputation tutorial.  
**Input file:**  
-- sample.vcf.gz  
**Outpt file:**  
-- sample_autosome.vcf.gz 
```bash
# Generate a comma-separated string of chromosome names to keep
chrs=$(paste -d ' ' <(echo {1..22}) | tr ' ' ',')

# Keep only those wanted chromosomes
bcftools view -t $chrs sample.vcf.gz -Oz -o sample_autosome.vcf.gz
```
> `paste` parameter:  
-d delimiter  
`bcftools view` parameters:  
-t targets  
^ exclusion prefix  
-Oz compressed output  

### 2.2. Check REF/ALT flips

Align the variant alleles to human reference genome to correct for any dataset-specific REF/ALT flips.  
Ensure that only biallelic sites are kept in the target data, as `bcftools norm` may introduce false multiallelic sites.  
Finally, replace the ID column with a 'SNP ID' in format CHROM_POS_REF_ALT ie. chromosome_position_<reference allele>_<alternative allele>.
 
**Input files:**  
-- sample_autosome.vcf.gz  
-- human_g1k_v37_decoy.fasta  
**Outpt file:**  
-- sample_fixmis.vcf.gz

```bash
# Align the alleles to the reference genome
bcftools norm -f hg38_v0_Homo_sapiens_assembly38.fasta -c ws <dataset>_chrfiltered.vcf.gz -Ou | \
# Keep only biallelic records
bcftools view -m 2 -M 2 -Oz -o <dataset>_refcorrected.vcf.gz

# Replace the ID column with a CHR_POS_REF_ALT 'SNP ID'
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' <dataset>_refcorrected.vcf.gz -Oz -o <dataset>_SNPID.vcf.gz
```
>`bcftools norm` parameters:  
-f reference genome  
-c what to do when incorrect/missing REF allele is encountered:  
ws warn (w) and set/fix (s) bad sites  
`bcftools view`  
-m minimum number of alleles listed in REF and ALT columns  
-M maximum number of alleles listed in REF and ALT columns  
-Ou uncompressed output  
-Oz compressed output  
`bcftools annotate` parameter:  
--set-id set ID column, % indicates a VCF field  


### 2.3. Remove duplicate samples 

Ensure that duplicate individuals do not exist in the chip data or between the chip and reference panel as this would compromise the imputation.
 
**Input files:**  
-- <dataset>_SNPID.vcf.gz  
-- panel_sample_IDs.txt  
**Outpt file:**   
-- <dataset>_noduplicate_samples.vcf.gz

```bash
# Copy the panel sample ID file with a different name to the working directory
cp /path/to/panel_sample_IDs.txt duplicate_sample_IDs.txt

# Generate a list of sample IDs from the chip data, keep only duplicates and append to the list of reference panel sample IDs
bcftools query -l <dataset>_SNPID.vcf.gz | uniq -d >> duplicate_sample_IDs.txt

# Remove the listed sample IDs from the chip data VCF
bcftools view -S ^duplicate_sample_IDs.txt --force-samples <dataset>_SNPID.vcf.gz -Oz -o <dataset>_noduplicate_samples.vcf.gz
```
>`bcftools query` parameters:  
-l list of sample IDs  
`uniq` parameter:  
-d only print duplicate lines  
`bcftools view` parameters:  
-S file of sample IDs to include  
^ exclusion prefix  
--force-samples only warn about unknown subset of samples  
-Oz compressed output  


### 2.4. Remove duplicate variants

Ensure that there are no duplicate variants (these might appear in some target genotype datasets).  
Duplicate variants will need to be removed before imputation.
 
**Input file:**  
-- <dataset>_noduplicate_samples.vcf.gz
**Outpt files:**  
-- <dataset>_duplicate_variants.txt  
-- <dataset>_noduplicate_variants.vcf.gz  
 
```bash
# Store a list of duplicate positions
bcftools query -f '%ID\n' <dataset>_noduplicate_samples.vcf.gz | uniq -d > <dataset>_duplicate_variants.txt

# Check whether the file contains any variants
if [ -s <dataset>_duplicate_variants.txt ]; then
   # Then remove the duplicate variants
   bcftools view -e ID=@<dataset>_duplicate_variants.txt <dataset>_noduplicate_samples.vcf.gz -Oz -o <dataset>_noduplicate_variants.vcf.gz
else
  # If the file is empty i.e. no duplicate variants are present, only rename the file to be compatible with the next step
  mv <dataset>_noduplicate_samples.vcf.gz <dataset>_noduplicate_variants.vcf.gz
fi
```
`bcftools query` parameters:  
-f query fields  
% identicate the field  
`uniq` parameter:  
-d only print duplicate lines  
`bcftools view` parameters:  
-e exclusion with expression  
`expression`:  
ID=@file IDs included in the file  
-Oz compressed output  