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
  make_option(c("-t", "--taxonomy"), action="store", default=NA, type='character',
              help=str_c("Taxonomy table with OTU IDs and taxonomic classification",
                         "(accepted extensions for this and other tabular data: .csv, .tab, .tsv)")),
  make_option(c("-o", "--otu"), action="store", default="OTU", type='character',
              help="'OTU' column name (default: OTU)"),
  make_option(c("-x", "--tax-columns"), action="store", default=NA, type='character',
              help="Comma-separated list of columns to retain in taxonomy table"),
  make_option(c("-O", "--otu-table"), action="store", default=NA, type='character',
              help="OTU table file (first column must contain OTU IDs)"),
  make_option(c("-f", "--fasta"), action="store", default=NA, type='character',
              help="OTU sequences file in fasta format (sequence IDs must match OTU names)"),
  make_option(c("-m", "--metadata"), action="store", default=NA, type='character',
              help="Metadata table (first column must contain sample IDs, which must match OTU table column names)"),          
  make_option(c("-c", "--cores"), action="store", default=NA, type='integer',
              help="Number of CPU cores"),
  make_option(c("-p", "--phyloseq"), action="store", default=NA, type='character',
              help="Output filename (.rds)"),
  make_option(c("-Z", "--optimize"), action="store_true", default=FALSE, type='logical',
              help="Optimize phylogenetic tree creation (will add considerable processing time)"),
  make_option(c("-T", "--no-tree"), action="store_true", default=FALSE, type='logical',
              help="Do not include a phylogenetic tree"))


#........................................................................
# Setup
#........................................................................

# parse command-line args
opt = parse_args2(OptionParser(option_list=option_list,formatter=nice_formatter))

taxonomy_file      <- opt$options$taxonomy
otu_file         <- opt$options$otu_table
otu_col       <- opt$options$otu
fasta_file    <- opt$options$fasta
metadata_file <- opt$options$metadata
cores         <- opt$options$cores
phyloseq_file <- opt$options$phyloseq
optimize      <- opt$options$optimize
no_tree       <- opt$options$no_tree
tax_columns   <- opt$options$tax_columns


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
if (!no_tree & !file_exists(fasta_file)) {
  stop(str_glue("Specified FASTA file '{fasta_file}' does not exist!"))
}

to_rename = c(otu = otu_col)
# load taxonomy table
taxonomy <- load_table(taxonomy_file,col_types=cols()) %>%
  dplyr::rename(all_of(to_rename)) %>%
  # select(where(is.character)) %>%
  select(otu,everything())

if (!is.na(tax_columns)) {
  to_keep <- tax_columns %>%
    str_split(",") %>%
    .[[1]]
  taxonomy <- taxonomy %>% select(otu,all_of(to_keep))
} else {
  taxonomy <- taxonomy %>% select(where(is.character))
}

# load otu table
otus <- load_table(otu_file,col_types=cols()) %>%
  dplyr::rename(otu=1)

# load sequences if they want a tree
if (!no_tree) {
  seqs <- read.fasta(fasta_file, as.string = TRUE, forceDNAtolower = FALSE) %>%
    map_dfr(~{
      c("otu" = attr(.x,"name"),"sequence"=.x)
    })
  taxonomy <- taxonomy %>%
    left_join(seqs %>% select(otu,sequence),by="otu") 
} else {
  taxonomy <- taxonomy %>%
    mutate(sequence=NA_character_)
}

tax_cols <- taxonomy %>%
  colnames()

otu_cols <- otus %>%
  colnames()

# do inner join so taxonomy and otus have same order ,etc.
combined <- otus %>%
  inner_join(taxonomy,by="otu")

# recreate original tables
otus <- combined %>%
  select(all_of(otu_cols))
taxonomy <- combined %>%
  select(all_of(tax_cols))
  

#........................................................................
# Prepare metadata
#........................................................................
# assume that the first column is the sample ID
metadata <- load_table(metadata_file,col_types=cols()) %>%
  column_to_rownames(colnames(.)[1])

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

if (!no_tree) {
  suppressPackageStartupMessages(library(DECIPHER))
  suppressPackageStartupMessages(library(phangorn))
  suppressPackageStartupMessages(library(Biostrings))
  # Phylogenetic tree code based on code from
  # https://ucdavis-bioinformatics-training.github.io/2021-May-Microbial-Community-Analysis/data_reduction/02-dada2
  dna <- DNAStringSet(taxonomy %>% pull(sequence))
  names(dna) <- taxonomy %>% pull(otu)

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
TREE   <- phy_tree(tree$tree)

if (!no_tree) {
  physeq <- phyloseq(OTU, TAX, META, TREE)
} else {
  physeq <- phyloseq(OTU, TAX, META)
}

saveRDS(physeq, file = paste0(phyloseq_file))