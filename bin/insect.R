#!/usr/bin/env Rscript

library(insect)
library(readr)
library(stringr)
library(purrr)

args <- commandArgs(TRUE)

zf <- args[1]
cf <- args[2]
cores <- as.numeric(args[3])
thresh <- as.numeric(args[4])
offs <- as.numeric(args[5])
mincount <- as.numeric(args[6])
pingr <- as.numeric(args[7])
output <- args[8]


zotus <- readFASTA(zf)
classifier <- readRDS(cf)
classified <- classify(
  zotus,
  classifier,
  threshold = thresh,
  metadata = TRUE,
  offset = offs,
  mincount = mincount,
  cores = cores,
  ping = pingr
)

write_csv(classified,output)

