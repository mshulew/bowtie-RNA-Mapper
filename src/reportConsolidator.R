#!/usr/bin/env Rscript

# combined multiple reports (report.tsv) into a single report
# arguments = parent output directory

library(readr)
library(dplyr)
library(stringr)

inputPath <- commandArgs(trailingOnly = TRUE)

dirs <- list.dirs(path=inputPath, full.names=TRUE, recursive=FALSE)

for(a in 1:length(dirs)){
  split_str <- dirs[a] %>%
    str_split("/") %>%
    .[[1]]
  sampleName <- split_str[length(split_str)]
  temp_df <- read_tsv(paste(dirs[a],"/report/report.tsv",sep=""), col_names=c("Field","Data"), comment = '##')
# replace Data cell with NA for rRNA Mapper
  temp_df$Data <- ifelse(grepl("# rRNA",temp_df$Field), NA, temp_df$Data)
  if(a == 1){
    consolidated_df <- temp_df %>%
      select(Field, Data) %>%
      rename(!!sampleName := Data)
  } else if(a > 1){
    consolidated_df <- consolidated_df %>%
      mutate(!!sampleName := temp_df$Data[match(Field,temp_df$Field)])
  }
}

# remove name of first column
colnames(consolidated_df)[1] <- NA

write_tsv(consolidated_df, 'consolidated_report.tsv', na="")
