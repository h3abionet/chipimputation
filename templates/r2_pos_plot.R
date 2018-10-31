#!/usr/bin Rscript
# 25.10.2018
# Daniel Schreyer
# R-squared - Position - Plot ####

# required packages
library(dplyr)
library(ggplot2)
library(optparse)
library(tidyr)

# takes input files as arguments
option_list <- list(
  make_option(c("-i", "--info"), action="store", default = NA, type='character',
              help = "Imputation .info file"),
  make_option(c("-t", "--target"), action="store", default = NA, type = 'character',
              help = "Target .frq file"),
  make_option(c("-o", "--output"), action="store", default = NA, type = 'character',
              help = "Output .png file"),
  make_option(c("-m", "--maf"), action="store", default = NA, type = 'character',
              help = "Minor allele frequency threshold in %")
)

args <- parse_args(OptionParser(option_list = option_list))

# Read in the allele frequency files
info <- read.table(args$i, sep = "\t", header = T)
frq <- read.table(args$t, sep = "\t", header = T)

# modify SNP ID of one file
frq <- frq %>% separate("SNP", c("CHR", "Position", "REF", "ALT"), "_")
frq <- frq %>% mutate(SNP = paste(frq$CHR, frq$POS, frq$REF, frq$ALT, sep = ":")) %>% select(c("CHR","POS","SNP"))

# merge tables
full <- merge(frq, dplyr::select(info, c("SNP", "Rsq", "Genotyped", "MAF", "ALT_Frq")), by = "SNP")

# change chromosome names
full$CHR <- paste("Chromosome",full$CHR, sep = " ")

# extract imputed and genotyped SNPs
Imputed <- filter(full, Genotyped == "Imputed")
Imputed <- filter(Imputed, Rsq != "-" | !is.na(Rsq))

genotyped <- filter(full, Genotyped == "Genotyped")

# display only  ~ 50,000 Imputed SNPs in the r2 - position plot
N <- ifelse(nrow(Imputed)> 50000, as.integer(nrow(Imputed)/50000, 1))
Imputed <- Imputed[seq(1, nrow(Imputed),N),]

# display only ~ 1,000 genotyped SNPs as ticks in the r2 - position plot
N2 <- ifelse(nrow(genotyped)>1000, as.integer(nrow(genotyped)/1000), 1)
genotyped <- genotyped[seq(1, nrow(genotyped),N2),]

# filter out MAF below a given value
AF_thresh <- ifelse(!is.na(args$m) & args$m > 0, args$m/100, 0)
Imputed <- filter(Imputed, MAF > AF_thresh)

# categorize MAF levels
Imputed <- mutate(Imputed, MAF2 = MAF)
Imputed$MAF2[Imputed$MAF2 < 0.005] <- "0.0001"
Imputed$MAF2[Imputed$MAF2 < 0.05 & Imputed$MAF2 >= 0.005] <- "0.005"
Imputed$MAF2[Imputed$MAF2 <= 0.5 & Imputed$MAF2 >= 0.05] <- "0.05"

# plot rsquared vs. SNP_position
r2_position_plot <- ggplot(data = Imputed, aes(x = POS, y = Rsq, color = MAF2)) +
  geom_point() +
  geom_hline(aes(yintercept = 0.3, color = "red"), show.legend = F) +
  labs(x = "Position [bp]", y = "Imputation accuracy (r-squared)") +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) + scale_x_continuous(labels = scales::comma) +
  facet_grid(~CHR) + theme_classic() +
  scale_colour_manual(name  ="MAF",
                      values=c("#999999", "#E69F00", "#56B4E9", "red"),
                      breaks=c("0.0001", "0.005", "0.05"),
                      labels=c("[0.0001,0.005)","[0.005, 0.05)","[0.05, 0.5]")) +
  geom_rug(data=genotyped, aes(x = POS), inherit.aes = F)

# save plot as a .png file
ggsave(file = args$o, r2_position_plot, width = 10, height = 5.5, units = "in")
