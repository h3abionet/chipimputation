#!/usr/bin Rscript
# 25.10.2018
# Daniel Schreyer
# Frequency Distribution between r2 and Info ####

# required packages
library(ggplot2)
library(dplyr)
library(optparse)
library(ggsci)

# takes input files as arguments
option_list <- list(
  make_option(c("-i", "--info"), action="store", default = NA, type = 'character',
              help = "Imputation .info files as a list"),
  make_option(c("-o", "--output"), action="store", default = NA, type = 'character',
              help = "Output .png file")
)
args <- parse_args(OptionParser(option_list = option_list))

# extract info files and reference panel name out of inputs
input <- args$i
inputs <- unlist(strsplit(input,","))

panels <- list()
for(panel in inputs){
  file <- unlist(strsplit(panel, "=="))[2]
  name <- unlist(strsplit(panel, "=="))[1]
  panels[paste0(name)] <- file
}

# read in .info files of each reference panel and merge them together in one table
i <- 1
for(file in panels){
  name <- names(panels)[i]
  panel <- read.table(as.character(file), sep = "\t", header = T)
  panel <- panel %>% mutate(R_Panel = paste0(name)) %>% select(c("SNP", "MAF","Rsq","Genotyped","R_Panel"))
  if(i > 1){
    full <- rbind(full, panel)
  }else{full <- panel}
  i <- i+1
}

# filter out non-imputed SNPs
Imputed <- filter(full, Genotyped == "Imputed")

# calculate mean Rsq and frequency of both reference panel
# MAF are rounded to 2 decimal places <- bin
Imputed <- Imputed %>% mutate( MAF = round(MAF, 2)) %>% group_by(R_Panel, MAF) %>% summarise(Rsq_mean = mean(Rsq), N = n())

# plot Frequency and Imputation Quality against MAF
p <- ggplot(Imputed, aes(x = MAF, color = R_Panel)) + 
  geom_point(aes(y = Rsq_mean*max(Imputed$N), shape = "Rsq"), size = 1) +
  geom_line(aes(y = Rsq_mean*max(Imputed$N))) +
  geom_point(aes(y = N, shape = "Frq"), size = 1) + 
  geom_line(aes(y = N))+
  scale_y_continuous(sec.axis = sec_axis(~./max(Imputed$N), name = "mean IQS per MAF bin",
                                         breaks = c(seq(0,1,0.2))), name = "SNPs count", 
                     labels = scales::comma) + labs(colour = "") + 
  geom_point(aes(x = MAF, y = Rsq_mean*max(Imputed$N), shape = "Rsq", color = R_Panel)
            , inherit.aes = F, size = 1) +
  geom_point(aes(y = N, shape = "Frq"), size = 1) +
  scale_color_npg(name = "Reference Panel") +
  scale_shape_manual(name = "", breaks = c("Rsq", "Frq"),
                    labels = c("mean IQS per MAF bin","SNPs count" ), values = c(18,4)) +
  theme_bw()
# save plot as .png file
ggsave(filename = args$o ,plot = p, width = 8, height = 5, units = "in")

