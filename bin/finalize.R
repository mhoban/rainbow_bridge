#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(fs))


load_table <- function(fn,col_types=cols(),...) {
  tabs <- c("tsv","tab","txt")
  commas <- c("csv")
  if (str_to_lower(path_ext(fn)) %in% commas){
    return(read_csv(fn,col_types=col_types,...))
  } else if (str_to_lower(path_ext(fn)) %in% tabs) {
    return(read_tsv(fn,col_types=col_types,...))
  } else {
    stop(str_glue("File '{fn}' must be either .csv or .tsv format"))
  }
  return(NULL)
}

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
  make_option(c("-R", "--remap"), action="store", default="", type="character", help="Taxonomy remapping table"),
  make_option(c("-F", "--filter"), action="store", default="", type="character", help="Taxonomy filtering table"),
  make_option(c("-i", "--insect"), action="store", default="", type="character", help="Taxonomy derived from insect classification"),
  make_option(c("-l", "--lca"), action="store", default="", type="character", help="Taxonomy derived from LCA classification"),
  make_option(c("-c", "--controls"), action="store", default="", type="character", help="File containing list of negative control samples"),
  make_option(c("-a", "--control-action"), action="store", default="remove", type="character", help="Action to perform on taxa found in negative control samples [default: remove]"),
  make_option(c("-d", "--decontam-method"), action="store", default="auto", type="character", help="Method to use when removing contamination with the `decontam` package"),
  make_option(c("-n", "--concentration"), action="store", default="", type="character", help="Tabular file containing DNA concentrations for each sample"),
  make_option(c("-T", "--control-threshold"), action="store", default=0.1, type="double", help="Minimum number of reads to allow from negative controls / threshold value when using `decontam` method"),
  make_option(c("-f", "--abundance-filter"), action="store_true", default=FALSE, type="logical", help="Filter zotus by minimum relative abundance"),
  make_option(c("-t", "--abundance-threshold"), action="store", default=0.0001, type="double", help="Minimum relative abundance for abundance filtering [default: 0.0001]"),
  make_option(c("-r", "--rarefy"), action="store_true", default=FALSE, type="logical", help="Rarefy read counts to minimum depth"),
  make_option(c("-m", "--rarefaction-method"), action="store", default="perm", type="character", help="Rarefaction method [default: perm]"),
  make_option(c("-u", "--permutations"), action="store", default=100, type="double", help="Number of permutations for permutational rarefaction [default: 100]"),
  make_option(c("-p", "--taxon-priority"), action="store", default="lca", type="character", help="Which method gets priority with differing taxonomies [default: lca]"),
  make_option(c("-N", "--min-reads"), action="store", default=1000, type="double", help="Minimum number of reads in a sample [default: 1000]"),
  make_option(c("-M", "--filter-min"), action="store_true", default=FALSE, type="logical", help="Filter out samples below minimum read threshold (step performed *before* rarefaction)"),
  make_option(c("-s", "--seed"), action="store", default=31337, type="double", help="Random number generator seed [default: 31337]"),
  make_option(c("-C", "--curated"), action="store", default="", type="character", help="LULU-curated zOTU table")
)

# parse command-line options
opt = parse_args(
  OptionParser(
    option_list=option_list,
    formatter=nice_formatter,
    prog="finalize.R",
    usage="%prog [options] <zotu_table>"
  ), 
  convert_hyphens_to_underscores = TRUE,
  positional_arguments = 1#,
  # args=debug_args
)


# get options
zotu_file        <- opt$args[1]
filter_file      <- opt$options$filter
remap_file       <- opt$options$remap
insect_file      <- opt$options$insect
lca_file         <- opt$options$lca
controls_file    <- opt$options$controls
control_action   <- opt$options$control_action
control_thresh   <- opt$options$control_threshold
decontam_method  <- opt$options$decontam_method
conc_file        <- opt$options$concentration
filter_abundance <- opt$options$abundance_filter
abundance_thresh <- opt$options$abundance_threshold
rarefy           <- opt$options$rarefy
rare_method      <- opt$options$rarefaction_method
rare_perm        <- opt$options$permutations
tax_priority     <- opt$options$taxon_priority
filter_min       <- opt$options$filter_min
min_reads        <- opt$options$min_reads
seed             <- opt$options$seed
curated          <- opt$options$curated

if (!file_exists(zotu_file)) {
  stop(str_glue("Missing/bad zotu table: {zotu_file}"))
}

if (! tax_priority %in% c("lca","insect")) {
  stop("--taxon-priority must be one of 'lca' or 'insect'")
}

if (! control_action %in% c("remove","subtract","decontam")) {
  stop("--control-action must be one of 'remove', 'subtract', or 'decontam'")
}

if (! rare_method %in% c("perm","phyloseq")) {
  stop("--rarefaction-method must be one of 'perm' or 'phyloseq'")
}

# first, merge taxonomies (if necessary)

# die if none of the files exist
if (all(!file_exists(c(lca_file,insect_file)))) {
  stop("Must provide at least one of insect and/or LCA taxonomies")
}

# load the taxonomies
lca <- NULL
insect <- NULL
if (file_exists(lca_file)) {
  lca <- load_table(lca_file)
}
if (file_exists(insect_file)) {
  insect <- load_table(insect_file)
}

cols <- c("domain","kingdom","phylum","class","order","family","genus","species")
sel <- c(str_c(cols,"_lca"),str_c(cols,"_insect"))

# both taxonomies were supplied
if (!is.null(lca) & !is.null(insect)) {
  taxonomy <- lca %>%
    full_join(insect,by="zotu",suffix=c("_lca","_insect")) %>%
    pivot_longer(all_of(sel),names_to = c("prefix","suffix"),names_pattern="(.*)_(insect|lca)") %>%
    pivot_wider(names_from = suffix, values_from = value) 
  if (tax_priority == "lca") {
    taxonomy <- taxonomy %>%
      mutate(taxon = coalesce(lca,insect),taxid = coalesce(taxid_lca,taxid_insect)) 
  } else {
    taxonomy <- taxonomy %>%
      mutate(taxon = coalesce(insect,lca),taxid = coalesce(taxid_insect,taxid_lca)) 
  }
  taxonomy <- taxonomy %>%
    select(-lca,-insect,-taxid_insect,-taxid_lca) %>%
    pivot_wider(names_from = prefix, values_from = taxon) %>%
    select(zotu,taxid,unique_hits,all_of(cols))
} else if (!is.null(lca)) {
  taxonomy <- lca %>%
    select(zotu,taxid,unique_hits,domain:species)
} else if (!is.null(insect)) {
  taxonomy <- insect %>%
    mutate(unique_hits=NA) %>%
    select(zotu,taxid,unique_hits,domain:species)
} else {
  stop("No taxonomy loaded, something went very wrong")
}

taxonomy_raw <- taxonomy

# remap taxonomy if requested
if (file_exists(remap_file)) {
  remap <- load_table(remap_file,col_select = c(level=1,val=2,newlevel=3,newval=4))

  # this next bit is very clever, but may be more than necessary
  # i.e., it could possibily be equally as clever, but shorter
  taxonomy <- taxonomy %>%
    pivot_longer(domain:species,names_to = "level") %>%
    left_join(remap,by=c("level","value" = "val")) %>%
    mutate(edit = case_when(
      !is.na(newlevel) & !is.na(newval) ~ TRUE,
      .default = FALSE
    )) %>%
    unite("level",c(level,newlevel),sep=";",na.rm=TRUE) %>%
    unite("value",c(value,newval),sep=";",na.rm=TRUE) %>%
    separate_longer_delim(c(level,value),delim = ";") %>%
    distinct(zotu,level,value,.keep_all = TRUE) %>%
    group_by(zotu,level) %>%
    filter(n() == 1 | edit) %>%
    ungroup() %>%
    select(-edit) %>%
    pivot_wider(names_from=level,values_from=value) %>%
    select(zotu,taxid,unique_hits,domain,kingdom,phylum,class,order,family,genus,species)  
}

# filter known contaminants, etc.
if (file_exists(filter_file)) {
  tax_filter <- load_table(filter_file,col_select = c(level=1,value=2,action=3))
  taxonomy <- taxonomy %>%
    pivot_longer(domain:species,names_to = "level") %>%
    left_join(tax_filter,by=c("level","value"))  %>%
    group_by(zotu) %>%
    mutate(keep = any(action == "retain"),filter=any(action == "filter")) %>%
    ungroup() %>%
    mutate(across(keep:filter,~replace_na(.x,FALSE))) %>%
    filter(keep) %>%
    filter(!filter) %>%
    select(-c(keep:filter),-action) %>%
    pivot_wider(names_from=level,values_from=value)
}


# join in zotu table, rename the first column to zotu
zotu_table <- load_table(zotu_file,col_select=c(zotu=1,everything())) %>%
  filter(zotu %in% taxonomy$zotu)

zotu_table_raw <- zotu_table

# filter negative controls
if (file_exists(controls_file)) {
  controls <- read_lines(controls_file)

  # sum reads in negative controls
  zotu_table <- zotu_table %>%
    mutate(ctrls = rowSums(pick(all_of(controls)),na.rm=T))
  zotu_table <- switch(
    control_action,
    remove = {
      zotu_table %>%
        filter(ctrls <= floor(control_thresh))
    },
    subtract = {
      zotu_table %>%
        mutate(
          across(-c(zotu,all_of(controls)),~.x-ctrls),
          across(-zotu,~replace(.x,which(.x < 0),0))
        )
    },
    decontam = {
      suppressPackageStartupMessages(library(decontam))

      zz <- zotu_table %>%
        select(-ctrls) %>%
        pivot_longer(-zotu,names_to="sample",values_to="reads") %>%
        pivot_wider(names_from="zotu",values_from="reads") %>%
        column_to_rownames("sample") %>%
        as.matrix()
      if (file_exists(conc_file)) {
        conc <- load_table(conc_file,col_select=c(sample=1,concentration=2)) %>%
          slice(match(rownames(zz),sample)) %>%
          pull(concentration)
      } else {
        conc <- NULL
      }
      negs <- rownames(zz) %in% controls
      contam <- isContaminant(zz,neg=negs,conc=conc,method=decontam_method,threshold=control_thresh)
      to_keep <- contam %>% as_tibble(rownames="zotu") %>% filter(!contaminant) %>% pull(zotu)
      zotu_table %>%
        filter(zotu %in% to_keep) 
    }
  ) %>%
    select(-ctrls,-all_of(controls)) %>%
    mutate(total = rowSums(pick(-zotu))) %>%
    filter(total > 0) %>%
    select(-total)
}

# filter by abundance
if (filter_abundance) {
  zt <- zotu_table %>%
    column_to_rownames("zotu") 
  totals <- colSums(zt)
  zotu_table <- zt %>%
    decostand("total",MARGIN=2) %>%
    mutate(across(everything(),~replace(.x,(.x / sum(.x)) < abundance_thresh,0) * totals[cur_column()])) %>%
    as_tibble(rownames="zotu") %>%
    mutate(reads = rowSums(pick(-zotu))) %>%
    filter(reads > 0) %>%
    select(-reads)
}


# filter by minimum abundance
if (filter_min) {
  zotu_table <- zotu_table %>%
    column_to_rownames("zotu") %>%
    select( where(~sum(.x) >= min_reads) ) %>%
    as_tibble(rownames="zotu")
}

# rarefy to minimum read depth
if (rarefy) {
  zotu_table <- switch(
    rare_method,
    perm = {
      suppressPackageStartupMessages(library(EcolUtils))

      set.seed(seed)
      zotu_table <- zotu_table %>%
        # column_to_rownames("zotu") %>%
        pivot_longer(-zotu,names_to="sample",values_to="reads") %>%
        pivot_wider(names_from="zotu",values_from="reads") %>%
        column_to_rownames("sample") %>%
        rrarefy.perm(n=rare_perm) %>%
        as_tibble(rownames="sample") %>%
        pivot_longer(-sample,names_to="zotu",values_to="reads") %>%
        pivot_wider(names_from="sample",values_from="reads")
    },
    phyloseq = {
      suppressPackageStartupMessages(library(phyloseq))
      zt <- zotu_table %>%
        column_to_rownames("zotu") %>%
        as("matrix")
      zotu_table <- zt %>%
        otu_table(taxa_are_rows=TRUE) %>%
        phyloseq() %>%
        rarefy_even_depth(rngseed=seed) %>%
        otu_table() %>%
        as.data.frame() %>%
        as_tibble(rownames="zotu")
    }
  )
}

write_tsv(zotu_table_raw,"zotu_table_raw.tsv",na="")
write_tsv(taxonomy_raw,"taxonomy_raw.tsv",na="")

final <- taxonomy %>%
  inner_join(zotu_table,by="zotu")
write_tsv(final,"zotu_table_final.tsv",na="")

if (file_exists(curated)) {
  # load curated zotu table and keep only zotus that ended up in the final shebang
  lulu_table <- load_table(curated,col_select=c(zotu=1,everything())) %>%
    filter(zotu %in% zotu_table$zotu) 
  final_curated <- taxonomy %>%
    inner_join(lulu_table,by="zotu")
  write_tsv(final_curated,"zotu_table_final_curated.tsv",na="")
}

