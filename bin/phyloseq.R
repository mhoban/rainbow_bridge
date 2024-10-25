#!/usr/bin/env Rscript
#===========================================================================
# Create phyloseq object (with phyloseq)
# based on here: https://joey711.github.io/phyloseq/import-data.html
# Seb Rauschert
# 22/07/2022
# Modified: Jessica Pearce
# 04/11/2022
# Modified: Adam Bennett
# 22/02/2023
# Modified: Mykle Hoban
# 11/07/2023 (That's in July), and subsequently multiple times thereafter
#===========================================================================

suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(fs))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(seqinr))

# help message formatter
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

load_table <- function(fn,...) {
  tabs <- c("tsv","tab","txt")
  commas <- c("csv")
  if (str_to_lower(path_ext(fn)) %in% commas){
    return(read_csv(fn,...))
  } else if (str_to_lower(path_ext(fn)) %in% tabs) {
    return(read_tsv(fn,...))
  } else {
    stop(str_glue("File '{fn}' must be either .csv or .tsv format"))
  }
  return(NULL)
}



# Define options for command line
option_list = list(
  make_option(c("-t", "--tree"), action="store_true", default=FALSE, type='logical',
              help="Generate a phylogenetic tree from OTU sequences"),
  make_option(c("-s", "--sequences"), action="store", default=NA, type='character',
              help="OTU sequences file in fasta format (sequence IDs must match OTU names)"),
  make_option(c("-c", "--cores"), action="store", default=parallel::detectCores(), type='integer',
              help="Number of CPU cores"),
  make_option(c("-o", "--out"), action="store", default=NA, type='character',
              help="Output filename (.rds)"),
  make_option(c("-z", "--optimize"), action="store_true", default=FALSE, type='logical',
              help="Optimize phylogenetic tree creation (will add considerable processing time)")
)


#........................................................................
# Setup
#........................................................................

# use debug arguments if we have 'em
if (exists('debug_args')) {
  opt_args <- debug_args
} else {
  opt_args <- commandArgs(TRUE)
}

# parse command-line args
opt <- parse_args(
  OptionParser(
    option_list=option_list,
    formatter=nice_formatter,
    prog="phyloseq.R",
    usage="%prog [options] <otu_table> <taxonomy> <metadata>"
  ), 
  convert_hyphens_to_underscores = TRUE,
  positional_arguments = 3,
  args = opt_args
)

otu_file      <- opt$args[1]
taxonomy_file <- opt$args[2]
metadata_file <- opt$args[3]
fasta_file    <- opt$options$sequences
cores         <- opt$options$cores
output_file   <- opt$options$out
optimize      <- opt$options$optimize
build_tree    <- opt$options$tree

# sanity checking: make sure files exist
if (!file_exists(taxonomy_file)) {
  stop(str_glue("Specified taxonomy table '{taxonomy_file}' does not exist!"))
}
if (!file_exists(otu_file)) {
  stop(str_glue("Specified otu table '{otu_file}' does not exist!"))
}
if (!file_exists(metadata_file)) {
  stop(str_glue("Specified metadata file '{metadata_file}' does not exist!"))
}
if (build_tree & !file_exists(fasta_file)) {
  stop(str_glue("Specified FASTA file '{fasta_file}' does not exist!"))
}

# load taxonomy table
taxonomy <- load_table(taxonomy_file,col_types=cols()) %>%
  rename(otu=1) 

# load otu table
otus <- load_table(otu_file,col_types=cols()) %>%
  dplyr::rename(otu=1)

sequences <- NULL
# load sequences if they want a tree
if (build_tree) {
  sequences <- read.fasta(fasta_file, as.string = TRUE, forceDNAtolower = FALSE) %>%
    map_dfr(~{
      c("otu" = attr(.x,"name"),"sequence"=.x)
    })
} 

tax_cols <- taxonomy %>%
  colnames()

otu_cols <- otus %>%
  colnames()

# do inner join so taxonomy and otus have same order ,etc.
combined <- otus %>%
  left_join(taxonomy,by="otu")

if (!is.null(sequences)) {
  # filter sequences by what's in the otu table/taxonomy and make sure
  # it's in the same order
  sequences <- sequences %>%
    filter(otu %in% combined$otu) %>%
    arrange(factor(otu,levels=combined$otu))
}

# recreate original tables
otus <- combined %>%
  select(all_of(otu_cols))
taxonomy <- combined %>%
  select(all_of(tax_cols))
  

#........................................................................
# Prepare metadata
# assume that the first column is the sample ID
# and make sure it has the sample sample IDs as the OTU table
#........................................................................
metadata <- load_table(metadata_file,col_types=cols()) %>%
  rename(sample=1) %>%
  filter(sample %in% (otus %>% select(-otu) %>% colnames())) %>%
  column_to_rownames("sample")

#........................................................................
# Prepare taxonomy data
#........................................................................
taxa <- taxonomy %>%
  column_to_rownames("otu") %>%
  as("matrix")


#........................................................................
# Prepare otu data
#........................................................................
otu <- otus %>%
  column_to_rownames("otu") %>%
  as("matrix")


#........................................................................
# Create phylogenetic tree
#........................................................................

if (build_tree) {
  suppressPackageStartupMessages(library(DECIPHER))
  suppressPackageStartupMessages(library(phangorn))
  suppressPackageStartupMessages(library(Biostrings))
  # Phylogenetic tree code based on code from
  # https://ucdavis-bioinformatics-training.github.io/2021-May-Microbial-Community-Analysis/data_reduction/02-dada2
  dna <- DNAStringSet(sequences %>% pull(sequence))
  names(dna) <- sequences %>% pull(otu)

  alignment <- AlignSeqs(dna, anchor=NA, processors=cores)

  phang_align <- phyDat(as(alignment, "matrix"), type="DNA")
  dm          <- dist.ml(phang_align)
  treeNJ      <- NJ(dm)

  fit    <- pml(treeNJ, data=phang_align)
  fitGTR <- update(fit, k=4, inv=0.2)

  if (optimize) {
    # Use 'optim.pml()' to optimize the model parameters
    fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
  }
  tree <- fitGTR
} else {
  tree <- NULL
}


#........................................................................
# Create phyloseq object
#........................................................................
OTU    <- otu_table(otu, taxa_are_rows = TRUE)
TAX    <- tax_table(taxa)
META   <- sample_data(metadata)

if (build_tree) {
  TREE   <- phy_tree(tree$tree)
  physeq <- phyloseq(OTU, TAX, META, TREE)
} else {
  physeq <- phyloseq(OTU, TAX, META)
}

saveRDS(physeq, output_file)
