#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(fs))
suppressPackageStartupMessages(library(stringr))

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
  make_option(c("-q", "--qcov"), action="store", default=NA, type='double', help="Minimum query coverage threshold"),
  make_option(c("-p", "--pid"), action="store", default=NA, type='double', help="Minimum percent match ID threshold"),
  make_option(c("-d", "--diff"), action="store", default=NA, type='double', help="Percent ID difference threshold for matching query coverage"),
  make_option(c("-u", "--keep-uncultured"), action="store_true", default=FALSE, type='logical', help="Retain uncultured/environmental/cloned, etc. sequences"),
  make_option(c("-i", "--intermediate"), action="store", default=NA, type='character', help="Store intermediate filtered BLAST results in specified file"),
  make_option(c("-s", "--semicolon"), action="store_true", default=FALSE, type='logical', help="Interpret taxids split by semicolon"),
  make_option(c("-m", "--merged"), action="store", default="", type="character", help="NCBI merged.dmp file"),
  make_option(c("-k", "--drop-blank"), action="store_true", default=TRUE, type="logical", help="Drop entries with completely blank taxonomic lineages"),
  make_option(c("-r", "--dropped"), action="store", default="dropped", type='character', help="Placeholder for dropped taxonomic levels (default: dropped)")
)

# parse command-line options
opt = parse_args(
  OptionParser(
    option_list=option_list,
    formatter=nice_formatter,
    prog="collapse_taxonomy.R",
    usage="%prog [options] <blast_result> <zotu_table> <lineage_dump> <output_table>"
  ), 
  convert_hyphens_to_underscores = TRUE,
  positional_arguments = 4
)

# check that passed files all exist and bail on failure
fe <- file_exists(opt$args[1:3])
if (any(!fe)) {
  bad <- fe[!fe]
  msg <- str_c(str_glue("Missing/bad filename: {names(bad)}"),collapse="\n")
  stop(msg)
}

# store filenames
blast_file <- opt$args[1]
zotu_table <- opt$args[2]
lineage_dump <- opt$args[3]
output_table <- opt$args[4]
qcov_thresh <- opt$options$qcov
pid_thresh <- opt$options$pid
diff_thresh <- opt$options$diff
filter_uncultured <- !opt$options$keep_uncultured
intermediate <- opt$options$intermediate
semicolon <- opt$options$semicolon
drop_blank <- opt$options$drop_blank
dropped <- opt$options$dropped
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

# read zotu table
zotus <- read_tsv(zotu_table, col_types = cols())

# perform initial filtering
filtered <- blast %>%
  # sort by descending pid
  arrange(desc(qcov),desc(pident)) %>%
  # keep only unique sequence matches since you can get a ton of matches to the same sequence
  # if it's something like a whole mito/genome 
  distinct(seqid,.keep_all = TRUE) %>%
  # add in count of unique blast hits
  add_count(zotu,name="unique_hits") %>%
  # keep only unique zotu/taxid combinations
  distinct(zotu,taxid,.keep_all = TRUE) %>%
  # filter by percentid and query coverage thresholds
  filter(pident >= pid_thresh & qcov >= qcov_thresh) 

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
filtered <- filtered %>%
  mutate(taxid=as.numeric(taxid)) %>%
  # now we create a sequential diff in percent id
  group_by(zotu) %>%
  # retain only results with the highest qcov (often this is just all 100%)
  filter(qcov == max(qcov)) %>%
  arrange(desc(qcov),desc(pident)) %>%
  mutate(diff=abs(c(0,diff(pident)))) %>%
  ungroup()

# successively reduce dataset by removing entries
# whose percent ID differs by the threshold
while (any(filtered$diff > diff_thresh)) {
  filtered <- filtered %>% 
    slice(-which(diff > diff_thresh)) %>%
    group_by(zotu) %>%
    arrange(desc(qcov),desc(pident)) %>%
    mutate(diff=abs(c(0,diff(pident)))) %>%
    ungroup()
}

if (!is.na(intermediate)) {
  write_tsv(filtered,intermediate,na="")
}

# NOTE: i don't know if there's any reason we can't just do it this way:
#       it avoids the loop and looks a lot cleaner
# filtered <- filtered %>%
  # group_by(zotu) %>%
  # filter(qcov == max(qcov)) %>%
  # arrange(desc(qcov),desc(pident)) %>%
  # mutate(diff=pident-pident[1]) %>%
  # filter(diff < diff_thresh)

  
# load the lineage dump, the underscores in col_types lets us skip columns
# because NCBI uses '\t|\t' as their delimiter, we'll have a bunch of columns
# that just contain a pipe
lineage <- read_tsv(lineage_dump,col_types = "i_c_c_c_c_c_c_c_c_c_",
                    col_names = c("taxid","species","fakespecies","genus","family","order","class","phylum","kingdom","domain")) %>%
  select(taxid,domain,kingdom,phylum,class,order,family,genus,species)

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
  left_join(lineage,by="taxid") 

# filter out any uncultured, environmental, etc.
if (filter_uncultured) {
  filtered <- filtered %>%
    filter(if_all(domain:species,~!replace_na(str_detect(.,"uncultured|environmental sample|clone|synthetic"),FALSE))) 
}

# drop hits with empty taxonomy (use in combination with merged to be sure)
if (drop_blank) {
  filtered <- filtered %>%
    filter(!if_all(domain:species,is.na))
}


# here we collapse taxonomic levels that differ across remaining blast results
filtered <- filtered %>%
  group_by(zotu) %>%
  summarise(across(domain:species,~if_else(n_distinct(.x) == 1,.x[1],dropped)),unique_hits=unique_hits[1]) %>%
  arrange(parse_number(zotu))

# join the zOTU table, using its first column as the zotu column
# use inner join because some zOTUs don't end up in the table
collapsed <- filtered %>%
  inner_join(zotus,by=set_names(colnames(zotus)[1],"zotu")) %>%
  select(domain:species,OTU=zotu,unique_hits,everything()) %>%
  # arrange(domain,kingdom,phylum,class,order,family,genus,species,parse_number(OTU))
  arrange(parse_number(OTU))

# save the output table
write_tsv(collapsed,output_table,na="")
