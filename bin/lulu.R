#!/usr/bin/env Rscript 
# This is a lulu script for post-clustering of Zotu table created by unoise3

suppressPackageStartupMessages(library(lulu))
suppressPackageStartupMessages(library(docopt))

# parse command-line options using docopt
doc <- "Usage:
  f.R [options] <zotu_table> <match_list> <curated_zotu_table> <curation_map> <rds>

Options: 
  -m <mr> --min-ratio <mr>         Minimim abundance ratio between a potential error and a potential parent to be identified as an error [default: 1]
  -t <mrt> --min-ratio-type <mrt>  Minimum ratio type. Accepted values: 'min', 'avg' [default: min]
  -a <mm> --min-match <mm>         Minimum threshold of sequence similarity for considering any OTU as an error of another [default: 84]
  -r <mrc> --min-rc <mrc>          Minimum relative co-ocurrence rate [default: 0.95]"
opt <- docopt(doc)

min_ratio <- opt[['min-ratio']]
min_ratio_type <- opt[['min-ratio-type']]
min_match <- opt[['min-match']]
min_rc <- opt[['min-rc']]
zotu_table <- opt[['zotu_table']]
match_list <- opt[['match_list']]
curated_zotu_table <- opt[['curated_zotu_table']]
curation_map <- opt[['curation_map']]
rds <- opt[['rds']]

# read zotu table and match list
zotus <- read.csv(zotu_table,sep='\t',header=TRUE,as.is=TRUE, row.names = 1)
matchlist <- read.table(match_list, header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)

# run lulu curation
curated_result <- lulu(
  zotus,
  matchlist,
  minimum_ratio_type = min_ratio_type,
  minimum_ratio = min_ratio,
  minimum_match = min_match,
  minimum_relative_cooccurence = min_rc
) 

lulu_table <- curated_result$curated_table

lulu_table <- data.frame(
  OTU = rownames(lulu_table),
  lulu_table,
  check.names = FALSE,
  row.names = NULL
)

# write curated zotu table
write.table(lulu_table,curated_zotu_table,sep="\t",row.names=FALSE,quote=FALSE)

# write curation map
write.table(curated_result$otu_map, curation_map, sep="\t",quote=FALSE)

# save R object
saveRDS(curated_result,rds)
