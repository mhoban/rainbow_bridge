#!/usr/bin/env Rscript

# TODO: add option to ignore singletons when collapsing

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(fs))
suppressPackageStartupMessages(library(stringr))

eq <- function(a,b) {
  an <- map_lgl(a,is.null)
  bn <- map_lgl(b,is.null)
  if (is.null(a) && is.null(b)) return(TRUE)
  if (is.null(a) && !is.null(b) || !is.null(a) && is.null(b)) return(FALSE)
  if (is.na(a) && is.na(b)) return(TRUE)
  if(is.na(a) && !is.na(b) || !is.na(a) && is.na(b)) return(FALSE)
  return(a==b)
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

na_opt <- function(opt, flag, val, parser, ...) {
  replace(val,val == "",NA) 
}

# set up option list
option_list <- list(
  make_option(c("-q", "--qcov"), action="store", default=NA, type='double', help="Minimum query coverage threshold"),
  make_option(c("-p", "--pid"), action="store", default=NA, type='double', help="Minimum percent match ID threshold"),
  make_option(c("-d", "--diff"), action="store", default=NA, type='double', help="Percent ID difference threshold for matching query coverage"),
  make_option(c("-f", "--filter-max-qcov"), action="store_true", default=FALSE, type='logical', help="Retain only records with the highest query coverage"),
  make_option(c("-t", "--taxon-filter"), action="callback", default=NA, type='character', help="Regex to filter taxa (e.g., uncultured/synthetic/environmental sequences)",callback=na_opt),
  make_option(c("-c", "--case-insensitive"), action="store_true", default=FALSE, type='logical', help="Perform case-insensitve taxon filtering"),
  make_option(c("-i", "--intermediate"), action="callback", default=NA, type='character', help="Store intermediate filtered BLAST results in specified file",callback=na_opt),
  make_option(c("-s", "--semicolon"), action="store_true", default=FALSE, type='logical', help="Interpret taxids split by semicolon"),
  make_option(c("-m", "--merged"), action="store", default="", type="character", help="NCBI merged.dmp file"),
  make_option(c("-k", "--drop-blank"), action="store_true", default=TRUE, type="logical", help="Drop entries with completely blank taxonomic lineages"),
  make_option(c("-r", "--dropped"), action="store", default="dropped", type='character', help="Placeholder for dropped taxonomic levels (default: dropped)"),
  make_option(c("-z","--zotu-table"),action="store",default=NA,type="character",help="Optional (tab-separated) OTU table to merge with results (first column must be OTU ID)")
)


# parse command-line options
opt <- parse_args(
  OptionParser(
    option_list=option_list,
    formatter=nice_formatter,
    prog="collapse_taxonomy.R",
    usage="%prog [options] <blast_result> <lineage_dump> <output_table>"
  ), 
  convert_hyphens_to_underscores = TRUE,
  positional_arguments = 4
)

# check that passed files all exist and bail on failure
fe <- file_exists(opt$args[1:2])
if (any(!fe)) {
  bad <- fe[!fe]
  msg <- str_c(str_glue("Missing/bad filename: {names(bad)}"),collapse="\n")
  stop(msg)
}

# store filenames
blast_file <- opt$args[1]
lineage_dump <- opt$args[2]
nodes_dump <- opt$args[3]
output_table <- opt$args[4]
qcov_thresh <- opt$options$qcov
pid_thresh <- opt$options$pid
diff_thresh <- opt$options$diff
taxon_filter <- opt$options$taxon_filter
ignore_case <- opt$options$case_insensitive
intermediate <- opt$options$intermediate
semicolon <- opt$options$semicolon
drop_blank <- opt$options$drop_blank
dropped <- opt$options$dropped
zotu_table_file <- opt$options$zotu_table

if (str_to_lower(dropped) == "na") {
  dropped <- NA_character_
}


# read blast results table
blast <- read_tsv(
  blast_file,
  col_names = c("zotu","seqid","taxid","species","commonname","kingdom","pident","length","qlen","slen","mismatch",
                "gapopen","gaps","qstart","wend","sstart","send","stitle","evalue","bitscore","qcov","qcovhsp"),
  col_types="ccccccnnnnnnnnnnncnnnn"
)

# CEB for multiple seqid matches, keep just the best one
best_score <- function(x,...) {
  x %>%
    group_by(...) %>%
    filter(evalue == min(evalue)) %>%  #CEB
    filter(pident == max(pident)) %>%  #CEB
    filter(bitscore == max(bitscore)) %>%  #CEB
    filter(mismatch == min(mismatch)) %>%  #CEB
    arrange(evalue,desc(pident),desc(bitscore),mismatch) %>%   #CEB
    slice(1) %>%
    ungroup()
}

# perform initial filtering
filtered <- blast %>%
  # sort by descending pid
  arrange(desc(qcov),desc(pident)) %>%
  # Collapse multiple seqid matches by best score combination
  # (per CEB's suggestion in #97)
  group_by(zotu,seqid) %>%
  summarise(
    across(taxid:kingdom,~.x[1]),
    pident=max(pident),
    across(length:slen,~.x[1]),
    mismatch=min(mismatch),
    across(gapopen:stitle,~.x[1]),
    evalue=min(evalue),
    bitscore=max(bitscore),
    qcov=max(qcov),
    qcovhsp=max(qcovhsp)
  ) %>%
  ungroup() %>%
  # add in count of unique blast hits
  add_count(zotu,name="unique_hits") %>%
  # keep only best zotu/taxid combinations
  best_score(zotu,taxid) %>%
  # filter by percent id and query coverage thresholds
  filter(pident >= pid_thresh & qcov >= qcov_thresh) %>%
  mutate(taxid = as.numeric(taxid))

if (semicolon) {
  # depending on how blast was run, we could potentially have multiple taxids,
  # so let's separate out rows based on the delimited ';'
  filtered <- filtered %>%
    separate_longer_delim(c(taxid,species,commonname),";") 
} else {
  # otherwise just get rid of them
  filtered <- filtered %>%
    filter(!str_detect(taxid,';'))
}

# load the lineage dump, the underscores in col_types lets us skip columns
# because NCBI uses '\t|\t' as their delimiter and we'd have a bunch of columns that just contain a pipe
lineage <- read_tsv(
  lineage_dump,
  col_types = "i_c_c_c_c_c_c_c_c_c_",
  col_names = c("taxid","taxon","species","genus","family","order","class","phylum","kingdom","domain")
) %>%
  select(taxid,domain,kingdom,phylum,class,order,family,genus,species,taxon)

# get new taxids for any merged taxa, if relevant
if (file_exists(opt$options$merged)) {
  merged <- read_tsv(opt$options$merged,col_types="i_i_",col_names=c("old_taxid","new_taxid"))
  filtered <- filtered %>%
    left_join(merged,by=c("taxid" = "old_taxid")) %>%
    mutate( taxid = coalesce(new_taxid,taxid)) %>%
  select(-new_taxid)
}

# join in ranked taxonomic lineage
filtered <- filtered %>%
  rename(blast_species=species,blast_kingdom=kingdom) %>%
  mutate(taxid=as.integer(taxid)) %>%
  left_join(lineage,by="taxid") %>%
  # this next line works because blast results always contain
  # species-level taxids
  mutate(species = coalesce(species,taxon)) 

# save intermediate file, if requested
if (!is.na(intermediate)) {
  write_tsv(filtered,intermediate,na="")
}

# now retain only assignments within the specified percent ID difference threshold
filtered <- filtered %>%
  # group by zotu
  group_by(zotu) %>%
  # conditionally retain only the highest query coverage within each zotu
  { if (opt$options$filter_max_qcov) filter(.,qcov == max(qcov)) else . } %>%
  # now calculate difference between each and the max pident within each zotu
  mutate(diff = abs(pident - max(pident))) %>%
  # discard anything with a difference over the threshold within each zotu
  filter(diff < diff_thresh) %>%
  ungroup() 

# filter taxa using supplied regex
if (!is.na(taxon_filter)) {
  filtered <- filtered %>%
    filter(if_all(domain:species,~!replace_na(str_detect(.,regex(taxon_filter,ignore_case = ignore_case)),FALSE))) 
}

# drop hits with empty taxonomy (use in combination with merged to be sure)
if (drop_blank) {
  filtered <- filtered %>%
    filter(!if_all(domain:species,is.na))
}

# here we collapse taxonomic levels that differ across remaining blast results
# keep taxids for species-level IDs, otherwise assign NA
filtered <- filtered %>%
  group_by(zotu) %>%
  summarise(
    across(domain:species,~if_else(n_distinct(.x) == 1,.x[1],dropped)),
    unique_hits=unique_hits[1],
    taxid = (\(tids) {
      if (n_distinct(tids) == 1) {
        return(unique(tids))
      } else {
        return(NA)
      }
    })(taxid)
  ) %>%
  arrange(parse_number(zotu))


# get taxid for last non-"dropped" taxonomic level
# skip this for now, because the join relationship goes wacky
# collapsed <- filtered %>%
#   mutate(
#     # first replace NAs with ... and dropped with NA (so coalesce works)
#     across(domain:species,~replace(.x,which(is.na(.x)),"...")),
#     across(domain:species,~replace(.x,which(.x == dropped),NA))
#   ) %>%
#   # now get the lowest non-NA taxonomic level
#   mutate(last_level = coalesce(species,genus,family,order,class,phylum,kingdom,domain)) %>%
#   # there is a problem when multiple names exist for something
#   # maybe use the actual level name in the join, if we can somehow
#   left_join(lineage,by=c("last_level" = "taxon"),relationship = "many-to-many",suffix=c("","_other")) %>%
#   mutate(
#     across(ends_with('_other'),~replace_na(.x,"")),
#     across(domain:species,~replace(.x,which(.x == "..."),"")),
#     across(domain:species,~replace_na(.x,dropped))
#   ) %>% add_count(zotu) %>% arrange(desc(n),zotu)

collapsed <- filtered

# merge zotu table, if desired
if (file_exists(zotu_table_file)) {
  zotu_table <- read_tsv(zotu_table_file,col_types=cols())
  # use inner join so we only get complete data
  collapsed <- collapsed %>%
    inner_join(zotu_table,by=setNames(colnames(zotu_table)[1],"zotu")) 
}

# save the output table
write_tsv(collapsed,output_table,na="")
