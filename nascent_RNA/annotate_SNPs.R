# Description:  In this script we will annotate the SNPs in the variants file with the dbSNP database

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
var_filename <- args[1]
dbsnp_filename <- args[2]

dbsnp <- as_tibble(readLines(dbsnp_filename)) %>%
  mutate(SNP = "SNP")
variants <- read.table(var_filename, header = TRUE, sep = "\t")

variants <- variants %>%
  mutate(id = paste(Region, Position, sep = ":"))

variants <- variants %>%
  left_join(dbsnp, by = c("id" = "value")) %>%
  filter(is.na(SNP)) %>%
  select(-SNP, -id)

outfile <- paste0(gsub(".txt", "", var_filename), "_noSNP.txt")

write.table(variants, file = outfile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")