#!/bin/bash

INPUT="./data"
OUTPUT="./output"

set -e

shapeit --input-bed ${INPUT}/gwas.bed ${INPUT}/gwas.bim ${INPUT}/gwas.fam \
        --input-map ${INPUT}/genetic_map.txt \
        --output-max ${OUTPUT}/gwas.phased.haps ${OUTPUT}/gwas.phased.sample


impute2 \
	-known_haps_g ${INPUT}/gwas_data_chr10_phased.haps \
	-h ${INPUT}/pilot1.jun2010.b36.CEU.chr10.snpfilt.haps \
	-l ${INPUT}/pilot1.jun2010.b36.CEU.chr10.snpfilt.legend \
	-m ${INPUT}/genetic_map_chr10_combined_b36.txt \
	-int 20000000 25000000 \ 
	-Ne 15000 \
	-buffer 250 \
	-o ${OUTPUT}/gwas_data_chr10_imputed.20-25Mb.gen
