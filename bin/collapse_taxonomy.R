#!/usr/bin/env Rscript

# TODO: add option to ignore singletons when collapsing

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(fs))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(R6))
suppressPackageStartupMessages(library(purrr))


# R6 class to get LCA of multiple NCBI taxids
# uses NCBI taxonomy dump, either nodes.dmp or taxidlineage.dmp
LCA <- R6Class(
  "LCA",
  private = list(
    merged = NULL,
    nodes = NULL,
    taxid_lineage = NULL,
    # get merged taxid
    get_merged = function(taxidlist) {
      if (!is.null(private$merged)) {
        taxidlist <- tibble(taxid=taxidlist) %>%
          left_join(private$merged,by=c("taxid" = "old_taxid")) %>%
          mutate(
            taxid = coalesce(new_taxid,taxid),
            name = if_else(!is.na(new_taxid),"merged",NA)
          ) %>%
          select(name,taxid) %>%
          tibble::deframe()
      }
      return(taxidlist)
    },
    valid = function(taxidlist) {
      if (!is.null(private$taxid_lineage)) {
        v <- taxidlist %in% private$taxid_lineage$taxid
      } else {
        v <- taxidlist %in% private$nodes$taxid
      }
      if (!is.null(private$merged)) {
        v <- v | v %in% private$merged$old_taxid
      }   
      return(v)
    }
  ),
  public = list(
    initialize = function(nodes = '',taxid_lineage='',merged = '') {
      # load nodes dump
      if (file_exists(nodes)) {
        private$nodes <- read_tsv(
          nodes,
          col_select = c(1,3,5),
          col_names=FALSE,
          progress=FALSE,
          show_col_types = FALSE
        ) %>%
          rename(taxid=1,parent=2,rank=3) %>%
          mutate(rank = replace(rank,which(rank == "superkingdom"),"domain"))
      } 
      
      # load merged taxid table
      if (file_exists(merged)) {
        private$merged <- read_tsv(
          merged,
          col_types="i_i_",
          col_names=c("old_taxid","new_taxid"),
          progress=FALSE,
          show_col_types = FALSE
        )
      }
      
      # load taxid lineage dump
      if (file_exists(taxid_lineage)) {
        private$taxid_lineage <- read_tsv(
          taxid_lineage,
          col_types="i_c_",
          col_names=c("taxid","lineage"),
          progress=FALSE,
          show_col_types = FALSE
        )
      }
      
      # throw error if neither nodes nor taxid lineage exist
      if (is.null(private$nodes) & is.null(private$taxid_lineage)) {
        stop("Must provide either nodes.dmp or taxidlineage.dmp")
      }
    },
    # get taxid lineage for a taxid
    get_lineage = function(tid) {
      if (!is.null(private$taxid_lineage)) {
        ll <- private$taxid_lineage %>%
          filter(taxid == tid)
        if (nrow(ll) > 0) {
          tids <- tibble(taxid=c(tid,rev(as.numeric(str_split_1(ll$lineage," "))))) %>%
            { 
              if (!is.null(private$nodes)) {
                left_join(.,private$nodes,by="taxid") %>%
                  select(rank,taxid) 
              } else .
            } 
          return(tids %>% tibble::deframe())
        } else {
          return(ll)
        }
      } else {
        node <- private$nodes %>%
          filter(taxid == tid[length(tid)])
        parent <- node %>%
          pull(parent)
        rank <- node %>%
          pull(rank)
        nn <- names(tid)
        if (length(nn) > 0) {
          nn[length(nn)] <- rank
        } else {
          nn <- rank
        }
        if (length(parent) > 0 & parent != 1) {
          return(self$get_lineage(c(setNames(tid,nn),parent)))
        } else {
          return(setNames(tid,nn))
        }
      }
    },
    lca = function(taxidlist, remove_bad_taxids = TRUE, ranks = "any") {
      if (remove_bad_taxids) {
        taxidlist <- taxidlist[private$valid(taxidlist)]
      } 
      
      taxidlist <- private$get_merged(taxidlist)  
      
      lineages <- taxidlist %>%
        map(\(tid) self$get_lineage(tid) )
      
      lcas <- lineages %>%
        reduce(intersect) 
      names(lcas) <- names(lineages[[1]])[which(lineages[[1]] %in% lcas)]
      
      if (all(ranks == "any") | all(is.null(names(lcas))) | all(names(lcas) == "")) {
        ranks <- 1
      } else {
        ranks <- ranks[na.omit(match(names(lcas),ranks))]
      }
      return(na.omit(lcas[ranks])[1])
    }
  )
)

# ordered hierarchy of taxa used by NCBI (and others)
# we should probably assemble this programatically, but they don't make it easy
hierarchy <- c( 
  "domain", "superkingdom", "supergroup", "kingdom", "subkingdom", "superphylum", "phylum", "division",
  "subphylum", "subdivision", "infraphylum", "superclass", "class", "subclass", "infraclass",
  "cohort", "subcohort", "superorder", "order", "suborder", "infraorder", "parvorder",
  "superfamily", "family", "subfamily", "tribe", "subtribe", "genus", "subgenus",
  "section", "subsection", "series", "subseries", "species group", "species subgroup",
  "species", "forma specialis", "subspecies", "varietas", "subvariety", "forma",
  "serogroup", "serotype", "strain", "isolate" 
)

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

# load a csv or tsv file
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

# validate NCBI lineage dump
check_ncbi_lineage <- function(f) {
  hdr <- read_lines(f,n_max = 1) %>%
    str_split_1("\t")
  if (length(hdr) == 20) {
    if (hdr[1] == "1" & hdr[3] == "root" & last(hdr) == "|") {
      return(TRUE)
    } 
  } 
  return(FALSE)
}


# set up option list
option_list <- list(
  make_option(c("-q", "--qcov"), action="store", default=NA, type='double', help="Minimum query coverage threshold"),
  make_option(c("-e", "--evalue"), action="store", default=NA, type='double', help="Maximum e-value threshold"),
  make_option(c("-p", "--pid"), action="store", default=NA, type='double', help="Minimum percent match ID threshold"),
  make_option(c("-d", "--diff"), action="store", default=NA, type='double', help="Percent ID difference threshold for matching query coverage"),
  make_option(c("-f", "--filter-max-qcov"), action="store_true", default=FALSE, type='logical', help="Retain only records with the highest query coverage"),
  make_option(c("-t", "--taxon-filter"), action="callback", default=NA, type='character', help="Regex to filter taxa (e.g., uncultured/synthetic/environmental sequences)",callback=na_opt),
  make_option(c("-c", "--case-insensitive"), action="store_true", default=FALSE, type='logical', help="Perform case-insensitve taxon filtering"),
  make_option(c("-i", "--intermediate"), action="callback", default=NA, type='character', help="Store intermediate filtered BLAST results in specified file",callback=na_opt),
  make_option(c("-s", "--semicolon"), action="store_true", default=FALSE, type='logical', help="Interpret taxids split by semicolon"),
  make_option(c("-m", "--merged"), action="store", default="", type="character", help="NCBI merged.dmp file"),
  make_option(c("-n", "--nodes"), action="store", default="", type="character", help="NCBI nodes.dmp file"),
  make_option(c("-a", "--taxid-lineage"), action="store", default="", type="character", help="NCBI taxidlineage.dmp file"),
  make_option(c("-k", "--drop-blank"), action="store_true", default=TRUE, type="logical", help="Drop entries with completely blank taxonomic lineages"),
  make_option(c("-r", "--dropped"), action="store", default="dropped", type='character', help="Placeholder for dropped taxonomic ranks (default: %default)"),
  make_option(c("-o", "--output"),action="store",default="lca_taxonomy.tsv",type="character",help="Output file (default: %default)")
)

# use debug arguments if we have 'em
if (exists('debug_args')) {
  opt_args <- debug_args
} else {
  opt_args <- commandArgs(TRUE)
}

# parse command-line options
opt <- parse_args(
  OptionParser(
    option_list=option_list,
    formatter=nice_formatter,
    prog="collapse_taxonomy.R",
    usage="%prog [options] <blast_result> <taxonomic_lineage>"
  ), 
  convert_hyphens_to_underscores = TRUE,
  positional_arguments = 2, 
  args = opt_args
)

# check that files in positional args all exist and bail on failure
fe <- file_exists(opt$args)
if (any(!fe)) {
  bad <- fe[!fe]
  msg <- str_c(str_glue("Missing/bad filename: {names(bad)}"),collapse="\n")
  stop(msg)
}

# get options
blast_file <- opt$args[1]
lineage_dump <- opt$args[2]
output_table <- opt$options$output
taxid_lineage_dump <- opt$options$taxid_lineage
nodes_dump <- opt$options$nodes
merged_dump <- opt$options$merged
qcov_thresh <- opt$options$qcov
eval_thresh <- opt$options$evalue
pid_thresh <- opt$options$pid
diff_thresh <- opt$options$diff
taxon_filter <- opt$options$taxon_filter
ignore_case <- opt$options$case_insensitive
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
  col_names = c("zotu","seqid","taxid","species","commonname","domain","pident","length","qlen","slen","mismatch",
                "gapopen","gaps","qstart","wend","sstart","send","stitle","evalue","bitscore","qcov","qcovhsp"),
  col_types="ccccccnnnnnnnnnnncnnnn",
  progress=FALSE,
  show_col_types = FALSE,
  na = "N/A"
)


# perform initial filtering
filtered <- blast %>%
  # sort by descending pid
  arrange(desc(qcov),desc(pident)) %>%
  # Collapse multiple seqid matches by best score combination
  # (per CEB's suggestion in #97)
  group_by(zotu,seqid) %>%
  summarise(
    across(taxid:domain,~first(.x)),
    pident=max(pident),
    across(length:slen,~first(.x)),
    mismatch=min(mismatch),
    across(gapopen:stitle,~first(.x)),
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
  { if (!is.na(pid_thresh)) filter(.,pident >= pid_thresh) else . } %>%
  { if (!is.na(qcov_thresh)) filter(.,qcov >= qcov_thresh) else . } %>%
  { if (!is.na(eval_thresh)) filter(.,evalue <= eval_thresh) else . } %>%
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


# lineage file is the NCBI dump
if (str_to_lower(path_file(lineage_dump)) == "rankedlineage.dmp") {
  if (check_ncbi_lineage(lineage_dump)) {
    # load the NCBI lineage dump, the underscores in lineage_cols lets us skip columns
    # because NCBI uses '\t|\t' as their delimiter and we'd have a bunch of columns that just contain a pipe
    # treat the 'taxon' column as the species column since blast always returns species-level taxids
    ncbi_ranks <- c("species","_","genus","family","order","class","phylum","kingdom","domain")
    lineage_cols <- str_c("i_", str_c( if_else(ncbi_ranks == "_","_","c"), collapse="_" ), "_")
    lineage_ranks <- ncbi_ranks[which(ncbi_ranks != "_")]
    lineage <- read_tsv(
      lineage_dump,
      col_types = lineage_cols,
      col_names = c("taxid",lineage_ranks),
      progress=FALSE,
      show_col_types = FALSE
    ) 
  } else {
    stop(str_glue("File {lineage_dump} is not a valid NCBI lineage dump"))
  }
} else {
  lineage <- load_table(lineage_dump,progress=FALSE,show_col_types=FALSE) %>%
    rename(taxid=1)
  if (!is.numeric(lineage$taxid)) {
    stop("First column of lineage table must be numeric taxon ID")
  }
  lineage_ranks <- colnames(lineage)[-1]
}


# get new taxids for any merged taxa, if relevant
if (file_exists(merged_dump)) {
  merged <- read_tsv(merged_dump,col_types="i_i_",col_names=c("old_taxid","new_taxid"),progress=FALSE, show_col_types=FALSE)
          
  filtered <- filtered %>%
    left_join(merged,by=c("taxid" = "old_taxid")) %>%
    mutate(taxid = coalesce(new_taxid,taxid)) %>%
    select(-new_taxid)
}

# join in ranked taxonomic lineage
filtered <- filtered %>%
  rename(blast_species=species,blast_domain=domain) %>%
  mutate(taxid=as.integer(taxid)) %>%
  left_join(lineage,by="taxid") 

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

# drop hits with empty taxonomy 
# (use in combination with merged.dmp to be sure you're not getting rid of real things)
if (drop_blank) {
  filtered <- filtered %>%
    filter(!if_all(all_of(lineage_ranks),is.na))
}

lca <- FALSE
# initialize LCA class if we're using it
# we do this if we have nodes.dmp or taxidlineage.dmp
if (file_exists(nodes_dump) | file_exists(taxid_lineage_dump)) {
  lca_getter <- LCA$new(nodes = nodes_dump, merged = merged_dump, taxid_lineage = taxid_lineage_dump)
  lca <- TRUE
}

# here we collapse taxonomic ranks that differ across remaining blast results
# keep taxids for species-level IDs
# otherwise (if we're using NCBI) assign taxid of LCA
collapsed <- filtered %>%
  group_by(zotu) %>%
  summarise(
    across(all_of(lineage_ranks),~ifelse(n_distinct(.x) == 1,first(.x),dropped)),
    unique_hits=unique_hits[1],
    taxid = (\(tids) {
      if (n_distinct(tids) == 1) {
        return(setNames(unique(tids),"species"))
      } else {
        if (lca) {
          # get taxid of LCA
          lt <- lca_getter$lca(tids,ranks = lineage_ranks)
          # make sure the resulting taxid has a name
          # either the name of its taxonomic rank or just a blank string
          return(setNames(lt,ifelse(all(is.null(names(lt))),"",names(lt))))
        } else {
          return(setNames(NA,""))
        }
      }
    })(taxid)
  ) %>%
  ungroup() %>%
  mutate(
    # assign taxonomic rank of resulting taxid using the name from above
    taxid_rank = case_when(
      length(names(taxid)) > 0 & names(taxid) != "" ~ names(taxid),
      !is.na(species) & species != dropped ~ "species",
      # length(names(taxid)) > 0 ~ names(taxid),
      .default = NA
    )
    ) %>%
  arrange(parse_number(zotu)) %>%
  # try to order the column hierarchically, matching the order of NCBI taxonomic hierarchy
  # if we have ranks not in the list, they'll end up at the end, but they'll still be there
  select(zotu,na.omit(match(hierarchy,colnames(.))),everything(),unique_hits,taxid,taxid_rank)


# save the collapsed output table
write_tsv(collapsed,output_table,na="")