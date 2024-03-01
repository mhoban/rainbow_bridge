#!/usr/bin/env Rscript 
# This is a lulu script for post-clustering of Zotu table created by unoise3

suppressPackageStartupMessages(library(lulu))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tibble))

nice_formatter <- function(object) {
  cat(object@usage, fill = TRUE)
  cat(object@description, fill = TRUE)
  cat("Options:", sep = "\n")

  options_list <- object@options
  for (ii in seq_along(options_list)) {
    option <- options_list[[ii]]
    cat("  ")
    if (!is.na(option@short_flag)) {
      cat(option@short_flag)
      if (optparse:::option_needs_argument(option)) {
        cat(" ", toupper(option@metavar), sep = "")
      }
      cat(", ")
    }
    if (!is.null(option@long_flag)) {
      cat(option@long_flag)
      if (optparse:::option_needs_argument(option)) {
        cat("=", toupper(option@metavar), sep = "")
      }
    }
    cat("\n    ")
    cat(sub("%default", optparse:::as_string(option@default), option@help))
    cat("\n\n")
  }
  cat(object@epilogue, fill = TRUE)
  return(invisible(NULL))
}

option_list = list(
  make_option(c("-m", "--min-ratio"), action="store", default=1, type='double', help="Minimum abundance ratio between a potential error and a potential parent to be identified as an error [default: 1]"),
  make_option(c("-t", "--min-ratio-type"), action="store", default="min", type='character', help="Minimum ratio type. Accepted values: 'min', 'avg' [default: min]"),
  make_option(c("-a", "--min-match"), action="store", default=84, type='double', help="Minimum threshold of sequence similarity for considering any OTU as an error of another [default: 84]"),
  make_option(c("-r", "--min-rc"), action="store", default=0.95, type='double', help="Minimum relative co-ocurrence rate [default: 0.95]")
)

# parse command-line options
opt = parse_args(
  OptionParser(
    option_list=option_list,
    formatter=nice_formatter,
    prog="lulu.R",
    usage="%prog [options] <zotu_table> <match_list> <curated_zotu_table> <curation_map> <rds>"
  ), 
  convert_hyphens_to_underscores = TRUE,
  positional_arguments = 5
)

# get command-line options
min_ratio          <- opt$options$min_ratio
min_ratio_type     <- opt$options$min_ratio_type
min_match          <- opt$options$min_match
min_rc             <- opt$options$min_rc
zotu_table         <- opt$args[1]
match_list         <- opt$args[2]
curated_zotu_table <- opt$args[3]
curation_map       <- opt$args[4]
rds                <- opt$args[5]

# read zotu table and match list
zotus <- read_tsv(zotu_table,col_types=cols()) %>%
  column_to_rownames(colnames(.)[1])  
matchlist <- read_tsv(match_list,col_names=c("child","parent","match"),col_types=cols())

# run lulu curation
curated_result <- lulu(
  zotus,
  matchlist,
  minimum_ratio_type = min_ratio_type,
  minimum_ratio = min_ratio,
  minimum_match = min_match,
  minimum_relative_cooccurence = min_rc
) 

lulu_table <- curated_result$curated_table %>% 
  as_tibble(rownames="zotu")
zotu_map <- curated_result$otu_map %>% 
  as_tibble(rownames="zotu")

# write curated zotu table
write_tsv(lulu_table,curated_zotu_table)

# write curation map
write_tsv(zotu_map,curation_map)

# save R object
saveRDS(curated_result,rds)
