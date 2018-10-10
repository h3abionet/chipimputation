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
singularity pull --name h3abionet-chipimputation-minimac4.simg shub://h3abionet/chipimputation:minimac4
```

Once transferred, start your singularity container:

```bash
singularity shell -B /data:/data  h3abionet-chipimputation-minimac4.simg
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
    wget -np ftp://share.sph.umich.edu/minimac3/G1K_P3_M3VCF_FILES_WITH_ESTIMATES.tar.gz
    ``` 
    
These phased reference panel VCF and M3VCF files are available in the shared folder `/data/refs/KGP` splitted in chromosomes.  
**NOTE**: The reference panel files should contain non-missing, phased genotypes!


## 2. Target/Chip data preparation

### 2.1. Get target data

For this tutorial, only a region on the genome will be used to reduce running time.
Extract chunk 6:8428721-8928720 for `sample.vcf.gz` with 5 samples only
```bash
# Index the vcf file
bcftools index --tbi -f sample.vcf.gz

# Generate a list of sample IDs present in the vcf file, one line per sample ID. Retaining only 5 samples
bcftools query -l sample.vcf.gz | head -n5 > panel_5sample_IDs.txt

# Subset the vcf only retaining the region of interest panel_5sample_IDs.txt
bcftools view --regions 6:8428721-8928720 --samples-file panel_5sample_IDs.txt --force-samples -m2 -M2 -v snps sample.vcf.gz -Oz -o sample5.vcf.gz
```
> `bcftools index` parameters:   
--tbi generate index for VCF file  
`bcftools view` parameters:  
--regions region to extract of format `chromosome:startPostion-endPosition  
-m2 minimum allele count   
-M2 maximum allele count   
-v type of variants to include 
--samples-list list of samples to include
--force-samples only warn about unknown subset of samples  
-Oz compressed output  

### 2.2. Check and remove invalid chromosomes

Only autosomal chromosomes are kept for this imputation tutorial.  
**Input file:**  
-- sample5.vcf.gz  
**Outpt file:**  
-- sample5_aut.vcf.gz 
```bash
# Generate a comma-separated string of chromosome names to keep
chrs=$(paste -d ' ' <(echo {1..22}) | tr ' ' ',')

# Keep only those wanted chromosomes
bcftools view -t $chrs sample5.vcf.gz -Oz -o sample5_aut.vcf.gz
```
> `paste` parameter:  
-d delimiter  
`bcftools view` parameters:  
-t targets  
^ exclusion prefix  
-Oz compressed output  

### 2.3. Check REF/ALT flips

Align the variant alleles to human reference genome to correct for any dataset-specific REF/ALT flips.  
Ensure that only biallelic sites are kept in the target data, as `bcftools norm` may introduce false multiallelic sites.  
Finally, replace the ID column with a 'SNP ID' in format CHROM_POS_REF_ALT ie. chromosome_position_<reference allele>_<alternative allele>.
 
**Input files:**  
-- sample5_chr6_8428721_8928720_aut.vcf.gz  
-- human_g1k_v37_decoy.fasta  
**Outpt file:**  
-- sample5_chr6_8428721_8928720_aut_fixmis.vcf.gz
-- sample5_chr6_8428721_8928720_aut_fixmis_SNPID.vcf.gz

```bash
# Align the alleles to the reference genome # Keep only biallelic records
bcftools +fixref sample5_aut.vcf.gz -Oz -o sample5_fix.vcf.gz -- -f /data/refs/common/human_g1k_v37_decoy.fasta -m flip -d 

# Replace the ID column with a CHR_POS_REF_ALT 'SNP ID'
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' sample5_fix.vcf.gz -Oz -o sample5_SNPID.vcf.gz
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


### 2.4. Remove duplicate samples 

Ensure that duplicate individuals do not exist in the chip data or between the chip and reference panel as this would compromise the imputation.  

First, we nned to generate a list of sample IDs from the reference panel, this can be done from any of the VCF files (here chr22) as in the example below (assuming that all chromosomes contain the same set of samples)

**Input files:**  
-- sample5_SNPID.vcf.gz  
-- /data/refs/KGP/vcf/1000GP_Phase3_chr6.vcf.gz  
**Output file:**   
-- panel_sample_IDs.txt
-- sample5_noduplicate_samples.vcf.gz

/*TODO*/ cp path
```bash
# Generate a list of sample IDs present in the reference panel, one line per sample ID
bcftools query -l /data/refs/KGP/vcf/1000GP_Phase3_chr6.vcf.gz > panel_sample_IDs.txt

# Copy the panel sample ID file with a different name
cp panel_sample_IDs.txt sample5_duplicate_sample_IDs.txt

# Generate a list of sample IDs from the chip data, keep only duplicates and append to the list of reference panel sample IDs
bcftools query -l sample5_SNPID.vcf.gz | uniq -d >> sample5_duplicate_sample_IDs.txt

# Remove the listed sample IDs from the chip data VCF
bcftools view -S ^sample5_duplicate_sample_IDs.txt --force-samples sample5_SNPID.vcf.gz -Oz -o sample5_noduplicate_samples.vcf.gz
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


### 2.5. Remove duplicate variants

Ensure that there are no duplicate variants (these might appear in some target genotype datasets).  
Duplicate variants will need to be removed before imputation.
 
**Input file:**  
-- sample5_noduplicate_samples.vcf.gz
**Outpt files:**  
-- sample5_duplicate_variants.txt  
-- sample5_noduplicate_variants.vcf.gz  
 
```bash
# Store a list of duplicate positions
bcftools query -f '%ID\n' sample5_noduplicate_samples.vcf.gz | uniq -d > sample5_noduplicate_samples.txt

# Check whether the file contains any variants
if [ -s sample5_duplicate_variants.txt ]; then
   # Then remove the duplicate variants
   bcftools view -e ID=@sample5_duplicate_variants.txt sample5_chr6_8428721_8928720_aut_fixmis_SNPID_noduplicate_samples.vcf.gz -Oz -o sample5_noduplicate_variants.vcf.gz
else
  # If the file is empty i.e. no duplicate variants are present, only rename the file to be compatible with the next step
  cp sample5_noduplicate_samples.vcf.gz sample5_noduplicate_variants.vcf.gz
fi
```
>`bcftools query` parameters:  
-f query fields  
% identicate the field  
`uniq` parameter:  
-d only print duplicate lines  
`bcftools view` parameters:  
-e exclusion with expression  
`expression`:  
ID=@file IDs included in the file  
-Oz compressed output  


### 2.6. Exclude rare variants

Generate a tab-delimited file of the target data allele frequencies, one variant per line, with columns CHR, SNP (in generated file header replaced with CHR_POS_REF_ALT), REF, ALT, AF (including the header line).  
First, exclude rare variants if not already removed in quality control steps, and re-calculate the allele frequency to correctly represent the current samples in the dataset.
 
**Input file:**  
-- sample5_noduplicate_variants.vcf.gz
 
**Output file:**  
-- sample5_AF.vcf.gz

```bash
# Remove low allele count variants if they are not already removed and re-calculate allele frequency
bcftools view -e 'INFO/AC<1 | (INFO/AN-INFO/AC)<1' sample5_noduplicate_variants.vcf.gz -Ou | \
bcftools +fill-tags  -Oz -o sample5_AF.vcf.gz -- -t AF
```
>`bcftools view` parameters:   
-e exclude based on expression  
-Ou uncompressed output  
`bcftools plugin` syntax and parameters:  
+fill-tags re-calculates/adds INFO field tags  
-Oz compressed output  
-- separator for plugin-specific parameters  
-t define the tags to be re-calculated/added. Alternatively, if all INFO field tags are wanted (see BCFtools documentation for complete list), remove the tag parameter: bcftools +fill-tags sample5_noduplicate_variants.vcf.gz -Oz -o sample5_noduplicate_variants_AF.vcf.gz


### 2.7. Generate target and reference panel allele frequencies

#### 2.7.1. Target data allele frequencies

Generate a frequency file for target data and reference panel allele frequency comparison

**Input file:**  
-- sample5_SNPID.vcf.gz
 
**Outpt file:**  
-- sample5_AF_chip.frq

```bash
# Re-calculate allele frequency
bcftools +fill-tags sample5_SNPID.vcf.gz -Oz -o sample5_AF.vcf.gz -- -t AF 

# First generate a tab-delimited header for the allele frequency file
echo -e 'CHR\tSNP\tREF\tALT\tAF' > sample5_AF_chip.frq

# Query the required fields from the VCF file and append to the allele frequency file
bcftools query -f '%CHROM\t%ID\t%REF\t%ALT\t%INFO/AF\n' sample5_AF.vcf.gz >> sample5_AF_chip.frq
```

`echo` parameter:  
-e enable interpretation of backslash escapes  
`bcftools query` parameter and syntax:  
-f format, where  
%field refers to a column in VCF file  
`expression` '%CHROM\t%ID\t%REF\t%ALT\t%INFO/AF\n' generates for each variant line tab-delimited format of: CHR, CHR_POS_REF_ALT, REF, ALT, AF


#### 2.7.2. Reference panel allele frequencies

Generate a tab-delimited file of the reference panel allele frequencies, one variant per line, with columns CHR, SNP (in generated file header replaced with CHR_POS_REF_ALT), REF, ALT, AF (including the header line); note last for loop in following Bash script.
 
Use the phased reference panel VCF files as input with the example command below and save the file as `panel.frq`

```bash
#Generate a tab-delimited header for the allele frequency file 
echo -e 'CHR\tSNP\tREF\tALT\tAF' > panel.frq  
```

```bash
# Check your reference panel VCF and if it does NOT contain AF in the INFO field, calculate it with `+fill-tags` plugin.   
#for chr in {1..23}; 
#do     
#    bcftools +fill-tags /data/refs/KGP/vcf/1000GP_Phase3_chr${chr}.vcf.gz -Oz -o /data/refs/KGP/vcf/1000GP_Phase3_chr${chr}_AF.vcf.gz -- -t AF 
#done  

# Query the required fields from the VCF file and append to the allele frequency file 
#for chr in {1..23}; 
#do     
#    bcftools query -f '%CHROM\t%CHROM\_%POS\_%REF\_%ALT\t%REF\t%ALT\t%INFO/AF\n' /data/refs/KGP/vcf/1000GP_Phase3_chr${chr}.vcf.gz >> panel.frq 
#done
```

For this tutorial, a preprocessed chromosome 6 reference panel VCF file with AF is located in `/data/imputation/KGP/1000GP_Phase3_chr6_AF.vcf.gz`.  
```bash
bcftools query -f '%CHROM\t%CHROM\_%POS\_%REF\_%ALT\t%REF\t%ALT\t%INFO/AF\n' /data/imputation/KGP/1000GP_Phase3_chr6_AF.vcf.gz >> panel.frq
```


### 2.8. Plot allele frequency differences between target and reference panel data and remove highly divergent variants

Compare the chip data and reference panel allele frequencies and generate an exclusion list of those variants where AF values differ more than 10 pp.
 
Download the R script given in 'FILE' link at the end of this step and save it as 'plot_AF.R'.
Set the script as executable for instance with 'chmod +x plot_AF.R' before running it for the first time.
 
Run the script as indicated in the example command by giving the two frequency files as input arguments.
 
**Input files:**  
-- sample5_AF_chip.frq  
-- panel.frq  
 
**Outpt files:**  
-- sample5_AFs.jpg  
-- sample5_AF_hist.jpg  
-- sample5_AF_exclude.txt 

```bash
# Copy and save the given 'plot_AF.R' file and run it with:
Rscript --no-save plot_AF.R sample5_AF_chip.frq panel.frq 
```

Save script to **plot_AF.R**
```R
#!/bin/env Rscript --no-save

# Genotype imputation protocol v3 - pre-imputation allele frequency plots
# Written by Kalle Pärn, Marita A. Isokallio, Paavo Häppölä, Javier Nunez Fontarnau and Priit Palta

# Required packages
library(data.table) # For fast fread()

# Input variables
args <- commandArgs(TRUE)
indata <- args[1]
paneldata <- args[2] 

# Read in the frequency files
chip <- fread(indata, header = T)
panel <- fread(paneldata, header = T)

# Generate a dataset tag
indataset <- sub("_chip.frq", "", indata)

# Take an intersection of the panel and chip data based on SNP column (in format CHR_POS_REF_ALT)
isec <- merge(panel, chip, by = "SNP")

# Check that AFs is within range of 10 pp in both datasets
af_ok <- abs(isec$AF.x - isec$AF.y) < 0.1

# Exclude those not within the AF range
exclude <- !af_ok

# Save the plot as jpg
jpeg(paste(indataset, "_AFs.jpg", sep=""))
# Plot first all and then excludable variants
plot(isec$AF.x, isec$AF.y, col=1, pch=20, main="Chip data AF vs. reference panel AF", xlab="Panel AF", ylab="Chip AF")
points(isec[exclude]$AF.x, isec[exclude]$AF.y, col=2, pch=20)
# Draw a legend
legend("topleft", legend=c("Concordant AF", "High AF difference"), col=c(1,2), pch=20, cex=0.9)
dev.off()

# Save the plot as jpg
jpeg(paste(indataset, "_AF_hist.jpg", sep = ""))
# Chip AF histogram for concordant AF variants
hist(isec[!exclude]$AF.y, breaks=100, main="Chip AF for concordant variants", xlab="Chip AF")
dev.off()

# Write out the exclusion list
write.table(isec[exclude]$SNP, paste(indataset, "_exclude.txt", sep=""), quote=F, row.names=F, col.names=F)
```

### 2.9. Exclude the variants with highly discordant variants.

Exclude the variants with highly discordant allele frequencies between the target dataset and the reference panel  as listed in the previous step.

**Input files:**   
-- sample5_AF_exclude.txt (output from `step 2.8`)  
-- sample5_AF.vcf.gz (output from `step 2.8`)  
 
**Outpt files:**  
-- sample5_AF_exclude_sorted.txt  
-- sample5_AF_for_phasing.vcf.gz  
 
```bash
# Sort the exclusion list
sort -V sample5_AF_exclude.txt > sample5_AF_exclude_sorted.txt

# Exclude the variants
bcftools view -e ID=@sample5_AF_exclude_sorted.txt sample5_AF.vcf.gz -Oz -o sample5_for_phasing.vcf.gz
tabix sample5_for_phasing.vcf.gz 
```

> `bcftools view` parameters:  
-e exclusion  
`expression`:  
ID=@file IDs in indicated file  
-Oz compressed output   


## 3. Chip/Target data phasing

Phasing of target data speeds up the genotype imputation step.  

**Input files:**  
-- sample5_for_phasing.vcf.gz  
-- /data/refs/KGP/genetic_map_hg19_withX.txt.gz  
-- /data/refs/KGP/vcf/1000GP_Phase3_chr6.vcf.gz  

**Outpt file:**  
-- sample5_for_imputation.vcf.gz  

```bash
# Run phasing for a specific region
eagle \
       --vcfTarget sample5_for_phasing.vcf.gz \
       --vcfRef=/data/refs/KGP/vcf/1000GP_Phase3_chr6.vcf.gz \
       --chrom 6 \
       --bpStart= 63854  \
       --bpEnd=20217978 \
       --bpFlanking=3000 \
       --geneticMapFile /data/refs/KGP/genetic_map_hg19_withX.txt.gz \
       --vcfOutFormat=z \
       --noImpMissing \
       --outPrefix sample5_for_imputation
```

Note: the 0|1, 1|1, et al represents phased genotypes, which are initially 0/1, 1/1 in the unphased data.
For example, for the first variant, 0 represents G, and 1 represents T.

> --vcf target data to phasing
--vcfRef reference panel to use for phasing
--chrom chromosome name
--geneticMapFile genetic map file provided by eagle
--vcfOutFormat output format to be compressed
--noImpMissing do not impute missing genotype
--outPrefix output prefix (to which vcf.gz will be added)

## 4. Genotype imputation

Run genotype imputation using minimac4.

**Input files:**  
-- sample5_for_imputation.vcf.gz  
-- /data/refs/KGP/m3vcf/1000GP_Phase3_chr6.m3vcf.gz
 
**Outpt file:**  
-- sample5_imputed.vcf.gz

```bash
 minimac4 \
        --refHaps /data/refs/KGP/m3vcf/1000GP_Phase3_chr6.m3vcf.gz \
        --haps sample5_for_imputation.vcf.gz \
        --format GT,DS \
        --allTypedSites \
        --minRatio 0.001 \
        --chr 6 --start 63854 --end 20217978 --window 3000 \
        --prefix sample5_imputed
```  
> `minimac4` parameters:
  --refHaps reference panel   
  --haps target data   
  --format GT,DS Produce genotype (GT) and (DS) allele dosage in the output vcf file  
  --allTypedSites Include genotyped variants in output  
  --minRatio 0.001 Ratio between variants in the target data vs the reference panel   
  --chr Chromosome to include in the analysis  
  --start Start position to consider in MB   
  --end End position to consider in MB  
  --window buffer size in MB  
  --prefix output file prefix  (prefix.dose.vcf.gz)
  
  
## 5. Post-imputation processing

Extract allele frequencies (AF) and INFO scores from the imputed VCF file.

Group the variants by their AF into either of the following categories,  
1: AF >= 5% or AF <= 95% (common variants)  
2: 0.5%  <= AF < 5% or 95% < AF <= 99.5% (low-frequency variants)  
3: AF < 0.5% or AF > 99.5% (rare variants)  
 
**Input file:**  
-- sample5_imputed.dose.vcf.gz   
 
**Outpt files:**  
-- sample5__varID_AF_INFO_GROUP.txt 
  
```bash
# Generate an allele frequency file for plotting for each chromosome

# Generate a header for the output file
echo -e 'CHR\tSNP\tREF\tALT\tAF\tINFO\tAF_GROUP' > sample5_varID_AF_INFO_GROUP.txt

# Query only the required fields and add allele frequency group (1, 2 or 3) as the ast column
bcftools query -f '%CHROM\t%CHROM\_%POS\_%REF\_%ALT\t%REF\t%ALT\t%INFO/AF\t%INFO/R2\t-\n' sample5_imputed.dose.vcf.gz |\
# $5 refers to AF values, $7 refers to AF group
awk -v OFS="\t" \
    '{if ($5>=0.05 && $5<=0.95) $7=1; \
       else if(($5>=0.005 && $5<0.05) || ($5<=0.995 && $5>0.95)) $7=2; \
       else $7=3} \
    { print $1, $2, $3, $4, $5, $6, $7 }' \
>> sample5_varID_AF_INFO_GROUP.txt

# Copy and save the given 'plot_INFO_and_AF_for_imputed_chrs.R' file and run it with:
Rscript --no-save plot_INFO_and_AF_for_imputed_chrs.R sample5 panel.frq

# Combine the plots per chromosome into a single pdf file
#convert $(ls <dataset>*postimputation_summary_plots.png | sort -V) sample5_postimputation_summary_plots.pdf
```
