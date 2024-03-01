#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(insect))
suppressPackageStartupMessages(library(fs))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))

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


# set up option list
option_list = list(
  make_option(c("-c", "--cores"), action="store", default=NA, type='double', help="Number of CPU cores to use"),
  make_option(
    c("-t", "--threshold"),
    action="store",
    default=FALSE,
    type='double',
    help="Minimum Akaike weight for the recursive classification procedure to continue toward the leaves of the tree"
  ),
  make_option(
    c("-o", "--offset"),
    action="store",
    default=NA,
    type='double',
    help="Log-odds score offset parameter governing whether the minimum score is met at each node"
  ),
  make_option(
    c("-n", "--min-count"),
    action="store",
    default=FALSE,
    type='double',
    help="Minimum number of training sequences belonging to a selected child node for the classification to progress"
  ),
  make_option(
    c("-p", "--ping"),
    action="store",
    default=NA,
    type="double",
    help="Value (between 0 and 1) indicating the minimum distance to the nearest neighbor for the the recursive classification algorithm to be skipped"
  ),
  make_option(c("-z","--zotu-table"),action="store",default=NA,type="character",help="Optional (tab-separated) OTU table to merge with results (first column must be OTU ID)"),
  make_option(c("-l","--lineage"),action="store",default=NA,type="character",help="NCBI ranked lineage file to fill in domain"),
  make_option(c("-m","--merged"),action="store",default=NA,type="character",help="NCBI merged taxid table")
)

# parse command-line options
opt = parse_args(
  OptionParser(
    option_list=option_list,
    formatter=nice_formatter,
    prog="insect.R",
    usage="%prog [options] <zotus> <model> <output>"
  ), 
  convert_hyphens_to_underscores = TRUE,
  positional_arguments = 3
)

# check that passed files all exist and bail on failure
fe <- file_exists(opt$args[1:2])
if (any(!fe)) {
  bad <- fe[!fe]
  msg <- str_c(str_glue("Missing/bad filename: {names(bad)}"),collapse="\n")
  stop(msg)
}

# pull out positional arguments
zotu_file <- opt$args[1]
model_file <- opt$args[2]
output_file <- opt$args[3]

# pull out options
cores <- opt$options$cores
thresh <- opt$options$threshold
offs <- opt$options$offset
mincount <- opt$options$min_count
pingr <- opt$options$ping
zotu_table_file <- opt$options$zotu_table
lineage_dump <- opt$options$lineage
merged_dump <- opt$options$merged

# read zotus
zotus <- readFASTA(zotu_file)
# read model
classifier <- readRDS(model_file)

# run the model
classified <- classify(
  zotus,
  classifier,
  threshold = thresh,
  metadata = TRUE,
  offset = offs,
  mincount = mincount,
  cores = cores,
  ping = pingr
) %>%
  mutate(taxID = as.numeric(taxID))

# get new taxids for any merged taxa, if relevant
if (file_exists(merged_dump)) {
  merged <- read_tsv(merged_dump,col_types="i_i_",col_names=c("old_taxid","new_taxid"))
  classified <- classified %>%
    left_join(merged,by=c("taxID" = "old_taxid")) %>%
    mutate( taxID = coalesce(new_taxid,taxID)) %>%
  select(-new_taxid)
}

# fill in domain column in classified output, if file is given
if (file_exists(lineage_dump)) {

  # load ranked lineage dump
  lineage <- read_tsv(
    lineage_dump,
    col_types = "i_c_c_c_c_c_c_c_c_c_",
    col_names = c("taxid","species","fakespecies","genus","family","order","class","phylum","kingdom","domain")
  ) %>%
    # just keep taxid and domain
    select(taxid,domain)
  
  # join lineage and make sure we're treaking Eukaryota as a domain
  classified <- classified %>%
    left_join(lineage,by=c("taxID" = "taxid")) %>%
    mutate(
      domain = case_when(
        taxon == "Eukaryota" ~ "Eukaryota",
        .default = domain
      )
    ) %>%
    # insert the domain column before the kingdom one
    relocate(domain,.before=kingdom)
}

# merge zotu table, if desired
if (file_exists(zotu_table_file)) {
  zotu_table <- read_tsv(zotu_table_file,col_types=cols())
  # use inner join so we only get complete data
  classified <- classified %>%
    inner_join(zotu_table,by=set_names(colnames(zotu_table)[1],"representative")) 
}

# save final output
write_tsv(classified,output_file)

