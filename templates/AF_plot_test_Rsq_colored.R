#!/usr/bin Rscript
# 25.10.2018
# Daniel Schreyer
# ref allele freq vs target allele freq Plot ####

# required packages
library(ggplot2)
library(dplyr)
library(optparse)
library(tidyr)

# takes Input files as arguments
option_list <- list(
  make_option(c("-i", "--info"), action="store", default=NA, type='character',
              help = "Imputation .info file"),
  make_option(c("-t", "--target"), action="store", default = NA, type = 'character',
              help = "Target .frg file"),
  make_option(c("-f", "--frq"), action="store", default = NA, type = 'character',
              help = "Reference Panel .frq file"),
  make_option(c("-r", "--rsq"), action="store", default = NA, type = 'character',
              help = "R-squared threshold"),
  make_option(c("-o", "--output"), action="store", default = NA, type = 'character',
              help = "Output .png file")
)
args <- parse_args(OptionParser(option_list = option_list))

# read in info and frequency file of imputed sample
info <- read.table(args$i, header = T, sep ="\t")
frq <- read.table(args$t, header = T, sep = "\t")

# modify SNP ID 
frq <- frq %>% separate("SNP", c("CHR", "Position", "REF", "ALT"), "_")
frq <- frq %>% mutate(SNP = paste(frq$CHR, frq$POS, frq$REF, frq$ALT, sep = ":")) %>% select(c("CHR","POS","SNP"))

# merge tables together and read in frequency file of reference panel
full <- merge(frq, dplyr::select(info, c("SNP","ALT_Frq", "Rsq", "Genotyped")), by = "SNP")
full <- merge(full, read.table(args$f , sep = "\t", header = T), by = c("CHR","POS"))

# filter rsquared SNPs above a given threshold and calculate the frequency difference
Rsq_thresh <- ifelse(!is.na(args$r), args$r, 0)
Imputed <- full%>% filter(Genotyped == "Imputed") %>% mutate(diff = abs(ALT_Frq-AF))

Imputed <- filter(Imputed, Rsq != "-" | !is.na(Rsq))

# filter every Nth SNP to extract 20000 SNPs --> Plot gets more clear
N <- ifelse(nrow(Imputed) > 20000,as.integer(nrow(Imputed)/20000),1)
Imputed <- Imputed[seq(1, nrow(Imputed),by = N),]

# plot the allele frequencies against each other
# SNP allele frequency difference greater than 0.15 are colored black
p <- ggplot(Imputed , aes(x = ALT_Frq, y = AF, color = Rsq)) +
  geom_point(size = 0.9) + theme_classic() + scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) + labs(x = "Ref Allele Frequency (Uploaded Samples)", 
                                                   y = "Ref Allele Frequency (Reference Panel)") +
  geom_text(color = "black", x = 0.2, y = 1.02, 
                                            label = paste("R-squared threshold:", Rsq_thresh, sep = " ")) +
  geom_text(color = "black", x = 0.2, y = 0.95, label = paste(nrow(Imputed),"SNPs", sep = " ")) +
  scale_color_gradient(low = "lightblue", high = "darkblue")

# save plot as a png file
ggsave(filename = args$o, plot = p, width = 7, height = 7, units = "in")


