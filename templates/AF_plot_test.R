#!/usr/bin Rscript
# 25.10.2018
# Daniel Schreyer
# ref allele freq vs target allele freq Plot ####

# required packages
library(ggplot2)
library(dplyr)
library(optparse)
library(tidyr)
library(data.table)

# takes Input files as arguments
option_list <- list(
  make_option(c("-i", "--info"), action="store", default="${info}", type='character',
              help = "Imputation .info file"),
  make_option(c("-t", "--target"), action="store", default ="${target}", type = 'character',
              help = "Target .frg file"),
  make_option(c("-f", "--frq"), action="store", default = "${ref_frq}", type = 'character',
              help = "Reference Panel .frq file"),
  make_option(c("-r", "--rsq"), action="store", default = "${Rsq_thresh}", type = 'character',
              help = "R-squared threshold"),
  make_option(c("-o", "--output"), action="store", default ="${plot_out}", type = 'character',
              help = "Output .png file")
)
args <- parse_args(OptionParser(option_list = option_list))

# read in info and frequency file of imputed sample 
# data.table function fread to select columns
info <- fread(as.character(args[1]), header = T, sep ="\t", select = c("SNP","ALT_Frq", "Rsq", "Genotyped"), stringsAsFactors = F)
frq <- fread(as.character(args[2]), header = T, sep = "\t", select = c("CHR","POS","SNP"), stringsAsFactors = F)

# modify SNP ID 
frq <- frq %>% separate("SNP", c("CHR", "Position", "REF", "ALT"), "_")
frq <- frq %>% mutate(SNP = paste(frq$CHR, frq$POS, frq$REF, frq$ALT, sep = ":"))

# merge tables together and read in frequency file of reference panel
full <- merge(frq, info, by = "SNP")
full <- merge(full, fread(as.character(args[3], select = c("CHR", "POS", "AF")), sep = "\t", header = T), by = c("CHR","POS"), stringsAsFactors = F)

# filter 
# rsquared SNPs above a given threshold and calculate the allele frequency difference
Rsq_thresh <- ifelse(!is.na(args[4]), as.numeric(args[4]), 0)
Imputed <- full %>% filter(Rsq != "-") %>% filter(Rsq >= as.numeric(Rsq_thresh) & Genotyped == "Imputed") %>%
  filter(!is.numeric(abs)) %>% mutate(diff = abs(ALT_Frq-AF))

# filter every Nth SNP to extract 20000 SNPs --> Plot gets more clear
N <- ifelse(nrow(Imputed) > 20000,as.integer(nrow(Imputed)/20000),1)
Imputed <- Imputed[seq(1, nrow(Imputed),by = N),]

# plot the allele frequencies against each other
# SNP allele frequency difference greater than 0.15 are colored black
p <- ggplot(filter(Imputed, diff <= 0.15) , aes(x = ALT_Frq, y = AF, color = diff)) +
  geom_point(size = 0.9) +theme_classic() + scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  scale_color_gradient(low = "darkblue", high = "lightblue") +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) + labs(x = "Ref Allele Frequency (Uploaded Samples)", 
                                                   y = "Ref Allele Frequency (Reference Panel)") +
  theme(legend.position="none") + geom_text(color = "black", x = 0.2, y = 1.02, 
                                            label = paste("R-squared threshold:", Rsq_thresh, sep = " ")) +
  geom_text(color = "black", x = 0.2, y = 0.95, label = paste(nrow(Imputed),"SNPs", sep = " ")) +
  geom_point(data = filter(Imputed, diff > 0.15), aes(x = ALT_Frq, y = AF),shape = 1, color = "black", size = 0.6)

# save plot as a png file
ggsave(filename = as.character(args[5]), plot = p, width = 7, height = 7, units = "in")