#!/usr/bin/env Rscript

# combined multiple gene_counts_rpkmtpm.tsv files and split read counts, rpkm and tpm into separate files
# arguments = parent output directory

library(readr)
library(dplyr)
library(tidyr)
library(stringr)

inputPath <- commandArgs(trailingOnly = TRUE)

countsfiles <- list.files(path=inputPath, full.names=TRUE, recursive=FALSE)

for(a in 1:length(countsfiles)){
  split_str <- countsfiles[a] %>%
    str_split("/") %>%
    .[[1]]
  
  split_str <- split_str[length(split_str)] %>%
    str_split("_") %>%
    .[[1]]
  
  sampleName <- split_str[1]
  
  temp_df <- read_tsv(countsfiles[a],col_types = cols(Geneid = col_character(),Chr = col_character())) %>%
    rename(Gene = Geneid)
  if(a == 1){
    counts_df <- temp_df %>%
      select(Gene, counts) %>%
      rename(!!sampleName := counts)
    rpkm_df <- temp_df %>%
      select(Gene, rpkm) %>%
      rename(!!sampleName := rpkm)
    tpm_df <-temp_df %>%
      select(Gene, tpm) %>%
      rename(!!sampleName := tpm)
  } else if(a > 1){
    counts_df <- counts_df %>%
      mutate(!!sampleName := temp_df$counts[match(Gene,temp_df$Gene)])
    rpkm_df <- rpkm_df %>%
      mutate(!!sampleName := temp_df$rpkm[match(Gene,temp_df$Gene)])
    tpm_df <- tpm_df %>%
      mutate(!!sampleName := temp_df$tpm[match(Gene,temp_df$Gene)]) 
  }
}

write_tsv(counts_df,'consolidated_counts.tsv')
write_tsv(rpkm_df,'consolidated_rpkm.tsv')
write_tsv(tpm_df,'consolidated_tpm.tsv')
