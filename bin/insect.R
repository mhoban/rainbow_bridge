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

settings <- c("zOTU file: {zf}","Classifier model: {cf}","cores: {cores}","threshold: {thresh}","offset: {offs}","mincount: {mincount}","ping: {pingr}")
settings <- map_chr(settings,str_glue)
write_lines(settings,"insect_settings.txt")


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

write_csv(classified,"insect_classified.csv")

