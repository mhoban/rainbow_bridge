#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// pull in the helper class
import helper
import colors

// save params to config file
def save_config(config_file) {
  // make sure the yaml is dumped in block format
  def opts = new org.yaml.snakeyaml.DumperOptions()
  opts.setDefaultFlowStyle(org.yaml.snakeyaml.DumperOptions.FlowStyle.BLOCK)
  opts.setPrettyFlow(true)
  def y = new org.yaml.snakeyaml.Yaml(opts)

  // dump the yaml file, converting durations to strings
  // because otherwise they get dumped weird
  new File(config_file).withWriter { w -> 
    y.dump(new LinkedHashMap(params).collectEntries{ [it.key.toString(), it.value instanceof nextflow.util.Duration ? it.value.toString() : it.value ]}, w) 
  }
}

// quick check if variable is numeric
def is_num(x) {
  return x instanceof Number
}

// quickly make a number
def num(x) {
  if (!is_num(x)) {
    return Float.parseFloat(x)
  } else {
    return x
  }
}

// sanity check to make sure command-line parameters are correct and valid
def check_params() {
  // show help message and bail
  if (params.help) {
    helper.usage(params)
    if (params.debug) {
      println("\n\n\n")
      println(params)
    }
    exit(0)
  }

  // check PCR primers
  if (params.fwdPrimer || params.reversePrimer) {
    // only allow primers or barcode
    if (params.barcode) {
      println(colors.red("Only one of ") + colors.bred("--barcode") + colors.red(" or ") + colors.bred("--fwd-primer/--reverse-primer") + colors.red(" may be passed"))
      exit(1)
    }
    // bail if either primer doesn't exist
    if ((!params.fwdPrimer) || (!params.reversePrimer)) {
      println(colors.red("Both forward and reverse primers are required"))
      exit(1)
    }
  }

  // check barcode file(s)
  if (params.barcode) {
    // only allow primers or barcode
    if (params.fwdPrimer || params.reversePrimer) {
      println(colors.red("Only one of ") + colors.bred("--barcode") + colors.red(" or ") + colors.bred("--fwd-primer/--reverse-primer") + colors.red(" may be passed"))
      exit(1)
    }
    // check barcode file(s) exist
    def f = file(params.barcode)
    f = helper.is_list(f) ? f : [f]
    if (!(f.size() && f.every{ it.exists() })) {
      println(colors.red("The specified barcode file(s) either don't exist or there was some problem"))
      exit(1)
    }
  } else {
    // bail if demultiplexed by barcode or pool and no barcode is given
    if (params.demultiplexedBy in ['pool','barcode']) {
      println(colors.red("A valid barcode file is required for this demultiplexing method"))
      exit(1)
    }
  }

  // bail if they have previously demultiplexed samples and they're trying to split them
  if (params.split && params.demultiplexedBy == "index") {
    println(colors.red("Parameters") + colors.bred(" --split ") + colors.red("and") + colors.bred(" index-based demultiplexing ") + colors.red("are mutually exclusive"))
    exit(1)
  }

  // make sure the right version of single,paired,demultiplexed is passed
  if (!helper.file_exists(params.sequences) && !params.standaloneTaxonomy && params.single == params.paired) {
    if (!params.single) {
      println(colors.red("One of either ") + colors.bred("--single") + colors.red(" or ") + colors.bred("--paired") + colors.red(" MUST be passed"))
    } else {
      println(colors.red("Only one of either ") + colors.bred("--single") + colors.red(" or ") + colors.bred("--paired") + colors.red(" may be passed"))
    }
    exit(1)
  }

  // validate phyloseq params
  if (params.phyloseq) {
    if (!helper.file_exists(params.metadata)) {
      println(colors.yellow("The specified phyloseql metadata file ('${params.metadata}') does not exist"))
      exit(1) 
    }

    switch(params.taxonomy) {
      case 'lca':
        if (!params.lca) {
          println(colors.yellow("You passed --phyloseq with 'lca' as the taxonomy option, but LCA has not been run."))
          println(colors.yellow("Did you forget the --lca option?"))
        }
        break
      case 'insect':
        if (!params.insect) {
          println(colors.yellow("You passed --phyloseq with 'insect' as the taxonomy option, but insect has not been run."))
          println(colors.yellow("Did you forget the --insect option?"))
        }
        break
      case 'combined':
        if (!params.lca && !params.insect) {
          println(colors.yellow("Note: phyloseq generation requires one of --insect, --lca, or a custom taxonomy table."))
        }
        break
      default:
        if (!helper.file_exists(params.taxonomy)) {
          println(colors.yellow("You passed --phyloseq with a user-supplied taxonomy table, but the file '${params.taxonomy}' does not exist"))
          /* exit(1) */
        }
        break
    }
  }

  // check to make sure standalone taxonomy will work
  if (params.standaloneTaxonomy) {
    if (!params.lca && !params.insect) {
      exit(1,colors.bred("--standalone-taxonomy") + colors.red(" requires ") + colors.bred("--insect") + colors.red(" and/or ") + colors.bred("--lca"))
    }

    if (params.lca) { 
      if (!helper.file_exists(params.blastFile)) {
        println(colors.red("The supplied blast result table \"${params.blastFile}\" does not exist"))
        exit(1)
      }
    }

    if (params.insect) {
      if (!helper.file_exists(params.insectSequences) && !helper.is_url(params)) {
        exit(1,colors.red("The supplied FASTA file \"${params.insectSequences}\" does not exist"))
      }
    }
  }

  // validate sample map
  if (params.sampleMap != "" && !helper.file_exists(params.sampleMap)) {
    println(colors.red("The supplied sample map file ${params.sampleMap} does not exist"))
    exit(1)
  }

  // check to make sure denoiser is a valid input
  if (!(params.denoiser in ['usearch','vsearch','dada2'])) {
    exit(1,colors.bred("--denoiser") + colors.red(" must be one of 'vsearch', 'usearch', or 'dada2'"))
  }

  // sanity check blast database
  if (params.blast) {
    // if --blast-taxdb is passed, check that it's a .tar.gz archive or true
    if (params.blastTaxdb && (params.blastTaxdb !== true && !(params.blastTaxdb =~ /(?i)\.tar\.gz$/) )) {
      println(colors.bred("--blast-taxdb") + colors.red(" must have no argument or be passed the path to a .tar.gz archive"))
      exit(1)
    }

    // exclude and include blast taxa are mutally exclusive, so let's
    // make sure only one of the options is passed (check --blastn-xxx options too)
    def pos = params.blastTaxa || params.containsKey('blastnTaxids')
    def neg = params.blastExcludeTaxa || params.containsKey('blastnNegative_taxids')
    if (pos && neg) {
      println( colors.red("Only one of ") + colors.bred("--blast-taxa") + colors.red(" or ") +
        colors.bred("--blast-exclude-taxa") + colors.red(" may be passed ") +
        colors.red("(this includes ") + colors.bred("--blastn-taxids") + colors.red(" and ") +
        colors.bred("--blastn-negative_taxids") + colors.red(")")
      )
      exit(1)
    }

    // get --blast param value(s) as a list with only unique elements
    def blasts = ([] + params.blast).unique()

    // make sure we've got at least one db
    if (!blasts.size()) {
      println(colors.red("You must pass at least one value to --blast"))
      exit(1)
    } else {
      // make sure all dbs exist
      blasts.each {
        if (!file("${it}.ndb").exists()) {
          println(colors.red("Could not find BLAST database '${it}'. Please provide the path to an existing blast database."))
          if (it =~ /~/) {
            println(colors.yellow("The BLAST database path '${it}' contains a tilde ('~') that was not expanded by the shell. Try entering an absolute path."))
          }
          exit(1)
        }
      }
    }
  } else {
    if (params.lca && !params.standaloneTaxonomy) {
      println(colors.red("--lca requires the --blast option."))
      exit(1)
    }
  }

  // check LCA options
  if (params.lca) {
    if (params.lcaDiff <= 0) {
      println(colors.red("--lca-diff argument must be a number greater than zero."))
      exit(1)
    }

    if (params.noTaxdump && !helper.file_exists(params.lcaLineage)) {
      println(
        colors.red("A custom lineage file (set with ") + colors.bred("--lca-lineage") +
        colors.red(") is required when ") + colors.bred("--no-taxdump") + colors.red(" is passed")
      )
      exit(1)
    }
  }

  // make sure insect parameter is valid: either a file or one of the pretrained models
  if (params.insect) {
    if (!helper.insect_classifiers.containsKey(params.insect.toLowerCase())) {
      if (!helper.file_exists(params.insect) && !helper.is_url(params.insect)) {
        println(
          colors.red("Value passed to ") + colors.bred("--insect") + 
          colors.red(" must be one of the supported builtins or an RDS file/URL") +
          colors.red(" containing a trained insect classifier model.")
        )
        println(colors.red("See rainbow_bridge.nf ") + colors.bred("--help") + colors.red(" for supported builtin models"))
        exit(1)
      }
    }
  }
}

/* DADA2 and related processes (these supercede some of the others below) */

// plot sequence read quality profiles
process dada_plot_quality_profiles {
  label 'r'
  label 'process_more_memory'

  publishDir "${params.outDir}/plots/quality"

  input:
    tuple val(key), path(reads)
  output:
    tuple val(key), path("${key}_quality_plot.pdf")
  
  script:
  """
  #!/usr/bin/env Rscript
  library(ggplot2)
  library(dada2)
  reads <- c(${reads.collect{ '"' + it + '"' }.join(",")})
  plotz <- plotQualityProfile(reads,n=${params.plotQualitiesN})
  ggsave(filename="${key}_quality_plot.pdf",plot=plotz,device=cairo_pdf,width=7,height=5,units="in")
  """
}

// trim primers from fwd/reverse reads
process trim_primers {
  label 'cutadapt'
  label 'process_low'

  publishDir "${params.preDir}/primers_trimmed"

  input:
    tuple val(key), path(reads), val(fwd), val(rev)
  output:
    tuple val(key), path("${key}_*_primer_trimmed.fastq.gz")

  script:
  // reverse complement the primers
  def bases = [
    'a': 'T', 't': 'A', 'u': 'A', 'g': 'C', 'c': 'G', 'y': 'R', 'r': 'Y', 's': 'S',
    'w': 'W', 'k': 'M', 'm': 'K', 'b': 'V', 'd': 'H', 'h': 'D', 'v': 'B', 'n': 'N',
    'A': 'T', 'T': 'A', 'U': 'A', 'G': 'C', 'C': 'G', 'Y': 'R', 'R': 'Y', 'S': 'S',
    'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N'
  ]
  def fwd_rc = fwd.reverse().collect{ bases[it] }.join('')
  def rev_rc = rev.reverse().collect{ bases[it] }.join('')
  """
  cutadapt \\
    --discard-untrimmed \\
    --no-indels \\
    -m ${params.minLen} \\
    -j ${task.cpus} \\
    -e ${params.primerMismatch} \\
    -a ${params.freePrimers ? "" : "^"}${fwd}...${rev_rc} \\
    -A ${params.freePrimers ? "" : "^"}${rev}...${fwd_rc} \\
    -o ${key}_R1_primer_trimmed.fastq.gz \\
    ${params.paired ? "-p ${key}_R2_primer_trimmed.fastq.gz" : "" } \\
    ${params.maxLen ? "-M "  + params.maxLen : '' } \\
    ${reads[0]} ${params.paired ? reads[1] : ""}
  """
}

// trim reads to minimum absolute length
process trim_length {
  label 'cutadapt'
  label 'process_low'

  publishDir "${params.preDir}/length_filtered"

  input:
    tuple val(key), path(reads)
  output:
    tuple val(key), path("${key}_*_length_trimmed.fastq.gz")

  script:
  """
  cutadapt \\
    -m ${params.minLen} \\
    -o ${key}_R1_length_trimmed.fastq.gz \\
    ${params.paired ? "-p ${key}_R2_length_trimmed.fastq.gz" : "" } \\
    ${reads[0]} ${params.paired ? reads[1] : ""}
  """
}

/* do dada2 filter and trim */
process dada_filter_trim {
  label 'r'
  label 'process_high'

  publishDir "${params.preDir}/filtered_trimmed", pattern: "*.fastq.gz"
  publishDir "${params.preDir}/filtered_trimmed", pattern: "*.tsv"

  input: 
    tuple val(samples), path(fwd), path(rev)

  output:
    tuple val(samples), path('*_R1_filtered_trimmed.fastq.gz'), path('*_R2_filtered_trimmed.fastq.gz'), path('fwd.rds'), path('rev.rds'), emit: result
    path('filter.rds'), emit: filter
    path('filter_report.tsv')
  
  script:
  if (params.paired) {
    """
    #!/usr/bin/env Rscript
    library(rlang)
    library(purrr)
    library(stringr)
    library(dada2)
    library(tibble)
    library(readr)
    library(dplyr)

    samples <- c(${samples.collect{"\"${it}\""}.join(",")})
    fwd <- c(${fwd.collect{"\"${it}\""}.join(",")})
    sample_map <- samples %>%
      set_names(fwd)  
    rev <- c(${rev.collect{"\"${it}\""}.join(",")})
    fwd_filtered <- map_chr(samples,\\(s) str_glue("{s}_R1_filtered_trimmed.fastq.gz")) %>%
      set_names(samples)
    rev_filtered <- map_chr(samples,\\(s) str_glue("{s}_R2_filtered_trimmed.fastq.gz")) %>%
      set_names(samples)

    filtered <- filterAndTrim(
      fwd = fwd,
      filt = fwd_filtered,
      rev = rev,
      filt.rev = rev_filtered,
      compress = TRUE,
      truncQ = ${params.dadaTruncQ},
      truncLen = c(${params.dadaTruncate ? params.dadaTruncate : params.dadaTruncF + "," + params.dadaTruncR}),
      trimLeft = ${params.dadaTrimLeft},
      trimRight = ${params.dadaTrimRight},
      maxLen = ${params.dadaMaxLen},
      minLen = ${params.dadaMinLen},
      maxN = ${params.dadaMaxN},
      maxEE = c(${params.dadaMaxEeF && params.dadaMaxEeR ? params.dadaMaxEeF + "," + params.dadaMaxEeR : params.dadaMaxEe + "," + params.dadaMaxEe}),
      rm.phix = TRUE,
      multithread = ${task.cpus},
      matchIDs = TRUE
    )

    kept <- filtered %>%
      as_tibble(rownames="file") %>%
      filter(reads.out > 0) %>%
      pull(file)
    fwd_filtered <- fwd_filtered[names(fwd_filtered) %in% sample_map[kept]]
    rev_filtered <- rev_filtered[names(rev_filtered) %in% sample_map[kept]]

    saveRDS(filtered,'filter.rds')

    filtered <- as_tibble(filtered,rownames="file")
    write_tsv(filtered,"filter_report.tsv")

    saveRDS(fwd_filtered,'fwd.rds')
    saveRDS(rev_filtered,'rev.rds')
    """
  } else if (params.single) {
    """
    #!/usr/bin/env Rscript
    library(rlang)
    library(purrr)
    library(stringr)
    library(dada2)
    library(tibble)
    library(readr)
    library(dplyr)

    samples <- c(${samples.collect{"\"${it}\""}.join(",")})
    fwd <- c(${fwd.collect{"\"${it}\""}.join(",")})
    sample_map <- samples %>%
      set_names(fwd)  
    fwd_filtered <- map_chr(samples,\\(s) str_glue("{s}_R1_filtered_trimmed.fastq.gz")) %>%
      set_names(samples)
    # make fake reverse reads
    rev_filtered <- map_chr(samples,\\(s) {
      fn <- str_glue("{s}_R2_filtered_trimmed.fastq.gz")
      system(paste("touch",fn))
      return(fn)
    }) %>%
      set_names(samples)

    filtered <- filterAndTrim(
      fwd = fwd,
      filt = fwd_filtered,
      compress = TRUE,
      truncQ = ${params.dadaTruncQ},
      truncLen = ${params.dadaTruncate},
      trimLeft = ${params.dadaTrimLeft},
      trimRight = ${params.dadaTrimRight},
      maxLen = ${params.dadaMaxLen},
      minLen = ${params.dadaMinLen},
      maxN = ${params.dadaMaxN},
      maxEE = ${params.dadaMaxEe},
      rm.phix = TRUE,
      multithread = ${task.cpus},
      matchIDs = TRUE
    )

    saveRDS(fwd_filtered,'fwd_unkept.rds')
    saveRDS(rev_filtered,'rev_unkept.rds')

    kept <- filtered %>%
      as_tibble(rownames="file") %>%
      filter(reads.out > 0) %>%
      pull(file)
    fwd_filtered <- fwd_filtered[names(fwd_filtered) %in% sample_map[kept]]
    rev_filtered <- rev_filtered[names(rev_filtered) %in% sample_map[kept]]

    saveRDS(filtered,'filter.rds')

    filtered <- as_tibble(filtered,rownames="file")
    write_tsv(filtered,"filter_report.tsv")

    saveRDS(fwd_filtered,'fwd.rds')
    saveRDS(rev_filtered,'rev.rds')
    """
  }
}

// learn error rate
process dada_learn_errors {
  label 'r'
  label 'process_medium'

  input:
    tuple path(reads), path(filtered), val(direction)
    // tuple val(samples), path(fwd), path(rev), path(image)
  output:
    tuple path(reads), path(filtered), path("error_${direction}.rds"), val(direction)
    // tuple val(samples), path(fwd), path(rev), path('errors.Rdata')

  script:
  """
  #!/usr/bin/env Rscript
  library(dada2)

  filtered <- readRDS("${filtered}")

  err <- learnErrors(filtered, multithread=${task.cpus})

  saveRDS(err,"error_${direction}.rds")
  """
}

process dada_plot_errors {
  label 'r'
  label 'process_more_memory'

  publishDir "${params.outDir}/plots/error"

  input:
    tuple path(reads), path(filtered), path(err), val(direction)
  output:
    path('*.pdf')

  script:
  """
  #!/usr/bin/env Rscript
  library(dada2)
  library(ggplot2)

  err <- readRDS("${err}")
  plotz <- plotErrors(err,nominalQ=TRUE)

  ggsave(filename="${direction}_error_plot.pdf",plot=plotz,device=cairo_pdf,width=6,height=6,units="in")
  """
}

// run the main dada2 sample inference algorithm
process dada_infer_samples {
  label 'r'
  label 'process_medium'

  input:
    tuple path(reads), path(filtered), path(error), val(direction)
  output:
    tuple path(reads), path(filtered), path(error), path("dada_${direction}.rds"), val(direction)

  script:
  """
  #!/usr/bin/env Rscript
  library(dada2)
  
  reads <- readRDS("${filtered}")
  err <- readRDS("${error}")

  dd <- dada(derep=reads,err=err,multithread=${task.cpus})

  saveRDS(dd,"dada_${direction}.rds")
  """
}

// merge forward and reverse ASVs
process dada_merge_reads {
  label 'r'
  label 'process_single'
  label 'process_more_memory'

  input:
    tuple path(reads), path(filtered), path(error), path(dada)
  output:
    path('merged.rds')

  script:
  """
  #!/usr/bin/env Rscript
  library(dada2)
  library(stringr)
  library(rlang)
  library(purrr)

  dir <- c('fwd','rev')
  dd <- map(set_names(dir,dir),\\(d) {
    list(
      filtered=readRDS(str_glue("{d}.rds")),
      error=readRDS(str_glue("error_{d}.rds")),
      dada=readRDS(str_glue("dada_{d}.rds"))
    )
  })

  merged <- mergePairs(
    dd\$fwd\$dada,dd\$fwd\$filtered,
    dd\$rev\$dada,dd\$rev\$filtered,
    verbose=TRUE
  )
  saveRDS(merged,'merged.rds')
  """
}

// generate ASV table and ASV fasta
process dada_make_asv_table {
  label 'r'
  label 'process_single'

  publishDir "${params.outDir}/asvs", pattern: "*.{tsv,fasta}"

  input:
    path(merged)
  output:
    path('asv_table.rds'), emit: asv
    path("asv_table.tsv")
    path("asvs.fasta")

  script:
  """
  #!/usr/bin/env Rscript
  library(dada2)

  merged <- readRDS('${merged}')
  asv_table <- makeSequenceTable(merged)
  sequences <- names(getUniques(asv_table))
  seqid <- paste0("ASV",seq(ncol(asv_table)))

  saveRDS(asv_table,'asv_table.rds')

  asv_table <- as.data.frame(t(asv_table))
  samples <- colnames(asv_table)
  asv_table\$ASV <- seqid
  rownames(asv_table) <- NULL
  asv_table <- asv_table[c('ASV',samples)]

  seqid <- paste0(">",seqid)
  fasta <- c(rbind(seqid,sequences))

  write.table(asv_table,"asv_table.tsv",sep="\\t",row.names=FALSE,quote=FALSE)
  writeLines(fasta,"asvs.fasta")
  """
}

// detect and remove chimeras
process dada_remove_chimeras {
  label 'r'
  label 'process_medium'

  publishDir "${params.outDir}/asvs", pattern: "*.{tsv,fasta}"

  input:
    path(asv_table)
  output:
    tuple path('asv_table_nochimeras.rds'), path(asv_table), emit: asv
    path('asvs_nochimeras.fasta'), emit: fasta
    path('asv_table_nochimeras.tsv'), emit: asv_table

  script:
  """
  #!/usr/bin/env Rscript
  library(dada2)

  asv_table <- readRDS("${asv_table}")
  asv_table_nochimeras <- removeBimeraDenovo(asv_table,method="${params.dadaChimeraMethod}",multithread=TRUE,verbose=TRUE)
  saveRDS(asv_table_nochimeras,'asv_table_nochimeras.rds')

  sequences <- names(getUniques(asv_table_nochimeras))
  seqid <- paste0("ASV",seq(ncol(asv_table_nochimeras)))

  asv_table_nochimeras <- as.data.frame(t(asv_table_nochimeras))
  samples <- colnames(asv_table_nochimeras)
  asv_table_nochimeras\$ASV <- seqid
  rownames(asv_table_nochimeras) <- NULL
  asv_table_nochimeras <- asv_table_nochimeras[c('ASV',samples)]

  seqid <- paste0(">",seqid)
  fasta <- c(rbind(seqid,sequences))

  write.table(asv_table_nochimeras,"asv_table_nochimeras.tsv",sep="\\t",row.names=FALSE,quote=FALSE)
  writeLines(fasta,"asvs_nochimeras.fasta")
  """
}

// output table tracking sequence lost across dada2 processes
process dada_track_reads {
  label 'r'
  label 'process_single'

  publishDir "${params.outDir}/asvs"

  input:
    tuple path(filter_table), path(filter), path(dada), path(merged), path(asv_table)
  output:
    path('dada_summary.tsv')
  
  script:
  if (params.paired) {
    """
    #!/usr/bin/env Rscript
    library(dada2)
    
    filter_table <- readRDS("${filter_table}")
    fwd_filtered <- readRDS('fwd.rds')
    rev_filtered <- readRDS('rev.rds')
    dada_fwd <- readRDS('dada_fwd.rds')
    dada_rev <- readRDS('dada_rev.rds')
    merged <- readRDS('merged.rds')
    asv_table_nochimeras <- readRDS('asv_table_nochimeras.rds')

    get_n <- function(x) sum(getUniques(x))
    sequence_table <- cbind(rownames(asv_table_nochimeras),filter_table, sapply(dada_fwd, get_n), sapply(dada_rev, get_n), sapply(merged, get_n), rowSums(asv_table_nochimeras))
    colnames(sequence_table) <- c("sample", "input", "filtered", "denoised_fwd", "denoised_rev", "merged", "chimeras_removed")

    write.table(sequence_table,"dada_summary.tsv",sep="\\t",row.names=FALSE,quote=FALSE)
    """
  } else {
    """
    #!/usr/bin/env Rscript
    library(dada2)
    
    filter_table <- readRDS("${filter_table}")
    fwd_filtered <- readRDS('fwd.rds')
    dada_fwd <- readRDS('dada_fwd.rds')
    asv_table_nochimeras <- readRDS('asv_table_nochimeras.rds')

    get_n <- function(x) sum(getUniques(x))
    sequence_table <- cbind(rownames(asv_table_nochimeras), filter_table, sapply(dada_fwd, get_n), rowSums(asv_table_nochimeras))
    colnames(sequence_table) <- c("sample", "input", "filtered", "denoised", "chimeras_removed")

    write.table(sequence_table,"dada_summary.tsv",sep="\\t",row.names=FALSE,quote=FALSE)
    """
  }
}

// trim and (where relevant) merge paired-end reads
process filter_merge {
  label 'process_medium'

  publishDir "${params.preDir}/trim_merge", mode: params.publishMode

  input:
    tuple val(key), path(reads)

  output:
    tuple val(key), path('*_trimmed_merged.fastq'), emit: result
    path 'settings.yml'

  script:
  if( params.single ) {
    // single end
    """
    echo "AdapterRemoval: \$(AdapterRemoval --version 2>&1 | awk '{print \$NF}')" >> settings.yml
    echo 'paired: false' >> settings.yml
    echo 'min-quality: ${params.minQuality}' >> settings.yml
    echo 'max-quality: ${params.maxQuality}' >> settings.yml
    echo 'mate-separator: ${params.mateSeparator}' >> settings.yml

    AdapterRemoval --threads ${task.cpus} --file1 ${reads[0]} \\
      --trimns --trimqualities \\
      --minquality ${params.minQuality} \\
      --qualitymax ${params.maxQuality} \\
      --mate-separator ${params.mateSeparator} \\
      --basename ${key}

    mv ${key}.truncated ${key}_trimmed_merged.fastq
    """
  } else if ( params.paired ) {
    // if reads are paired-end then merge
    """
    echo "AdapterRemoval: \$(AdapterRemoval --version 2>&1 | awk '{print \$NF}')" >> settings.yml
    echo 'paired: true' >> settings.yml
    echo 'min-quality: ${params.minQuality}' >> settings.yml
    echo 'max-quality: ${params.maxQuality}' >> settings.yml
    echo 'min-align-len: ${params.minAlignLen}' >> settings.yml
    echo 'mate-separator: ${params.mateSeparator}' >> settings.yml

    AdapterRemoval --threads ${task.cpus} --file1 ${reads[0]} --file2 ${reads[1]} \\
      --collapse --trimns --trimqualities \\
      --minquality $params.minQuality \\
      --qualitymax ${params.maxQuality} \\
      --minalignmentlength ${params.minAlignLen} \\
      --mate-separator ${params.mateSeparator} \\
      --basename ${key}

    mv ${key}.collapsed ${key}_trimmed_merged.fastq
    """
  }
}

process filter_ambiguous_indices {
  label 'obitools'
  label 'process_single'

  publishDir "${params.preDir}/index_filtered", mode: params.publishMode

  input:
    tuple val(key), path(reads)

  output:
    tuple val(key), path("*_valid_index.fastq")

  script:
  """
  obigrep --uppercase -D ':[ACGT]+\\+[ACGT]+\$' ${reads} > "${key}_valid_index.fastq"
  """
}

// replace I's with N's in barcode file primer sequences
process fix_barcodes {
  label 'shell'
  label 'process_single'

  input:
    path(barcodes)

  output:
    path("${barcodes.baseName}_fixed.${barcodes.extension}")

  script:
  """
  if [[ -e "${barcodes}" ]]; then
    fix_barcode.awk "${barcodes}" > "${barcodes.baseName}_fixed.${barcodes.extension}"
  else
    touch "${barcodes.baseName}_fixed.${barcodes.extension}"
  fi
  """
}

// split and properly modify barcode file if they're pooled
process split_barcodes {
  label 'shell'
  label 'process_single'

  input:
    path(barcodes)

  output:
    path('bc/*.tsv')

  script:
  """
  mkdir bc
  split_barcode.awk -v parent="${barcodes.baseName}" ${barcodes}
  """
}

// primer mismatch & sample assignment
// multiple different barcode files are possible
process ngsfilter {
  label 'obitools'
  label 'process_single'
  label 'process_more_memory'

  publishDir "${params.preDir}/ngsfilter", mode: params.publishMode

  input:
    tuple val(key), path(read), path(barcode)


  output:
    tuple val(key), path("*_annotated.fastq"), val("${barcode.baseName}"), emit: result
    path 'settings.yml'

  script:
  """
  echo "obitools: 1.2.13" >> settings.yml
  echo 'primer-mismatch: ${params.primerMismatch}' >> settings.yml

  ngsfilter --uppercase -t ${barcode} -e ${params.primerMismatch} -u "${key}_filter_orphans.fastq" ${read} > "${key}_${barcode.baseName}_annotated.fastq"
  """
}

// combine outputs from (possible) multiple barcode files, filter by length
process filter_length {
  label 'obitools'
  label 'process_single'

  publishDir "${params.preDir}/length_filtered", mode: params.publishMode

  input:
    tuple val(key), path(fastq), val(barcode)

  output:
    tuple val(key), path('*_length_filtered.fastq'), val(barcode), emit: result
    path 'settings.yml'

  script:
  """
  echo "obitools: 1.2.13" >> settings.yml
  echo 'min-len: ${params.minLen}' >> settings.yml

  obigrep --uppercase -l ${params.minLen} "${fastq}" > "${key}_length_filtered.fastq"
  """
}

// for non-demultiplexed runs, split the annotated reads file by samples
process split_samples {
  label 'obitools'
  label 'process_single'

  publishDir "${params.preDir}/split_samples", mode: params.publishMode

  input:
    tuple val(key), path('to_split'), val(barcode)

  output:
    path('*.fastq'), optional: true

  script:
  """
  obisplit --uppercase -t sample -u "${key}.orphans" to_split
  """
}

// relabel files for dereplication
process relabel {
  label 'denoiser'
  label 'process_low'

  publishDir "${params.preDir}/relabeled", mode: params.publishMode

  input:
    tuple val(key), path(fastq)
  output:
    path('*_relabeled.fasta'), optional: true, emit: result
    path 'settings.yml'


  script:
  // we have to convert everything to uppercase because obisplit --uppercase is broken
  // and usearch -otutab will treat lowercase sequences as masked
  // vsearch might as well, so we play it safe
  if (params.denoiser == "vsearch") {
    """
    echo "vsearch: \$(vsearch --version 2>&1| head -1 | awk  '{print \$2}' | sed 's/,\$//')" >> settings.yml
    echo 'denoiser: vsearch' >> settings.yml

    vsearch --threads ${task.cpus} --fastq_qmax ${params.maxQuality} --fastx_filter ${fastq} --relabel "${key}." --label_suffix ";sample=${key}" --fastaout - | \\
      awk '/^>/ {print;} !/^>/ {print(toupper(\$0))}' > "${key}_relabeled.fasta"
    """
  } else {
    """
    echo "usearch: \$(usearch | head -1 | awk '{print \$2}')" >> settings.yml
    echo 'denoiser: usearch' >> settings.yml

    # usearch doesn't allow output to stdout so we have to use an intermediate file
    usearch -fastq_filter ${fastq} -relabel "${key}." -fastaout tmp.fasta  -sample "${key}"
    awk '/^>/ {print;} !/^>/ {print(toupper(\$0))}' tmp.fasta > "${key}_relabeled.fasta"
    rm tmp.fasta
    """
  }
}

// concatenate all relabeled files. we only do this in a process
// instead of using collectFile so we can see that it's happening
// in the process list. also, the output from collectFile may not be cached and this will be
process merge_relabeled {
  label 'shell'
  label 'process_single'

  publishDir "${params.preDir}/merged"

  input:
    path('input-fastq?????.fastq')

  output:
    path("${params.project}_relabeled_merged.fasta")

  script:
  """
  cat input-*.fastq > "${params.project}_relabeled_merged.fasta"
  """
}

// dereplicate to unique sequences
process dereplicate {
  label 'denoiser'
  label 'process_full'

  publishDir "${params.outDir}/zotus", mode: params.publishMode

  input:
    tuple val(id), path(relabeled_merged)

  output:
    tuple val(id), path("${id}_unique.fasta"), emit: result
    path 'settings.yml'

  script:
  if (params.denoiser == "vsearch") {
    """
    echo "vsearch: \$(vsearch --version 2>&1| head -1 | awk  '{print \$2}' | sed 's/,\$//')" >> settings.yml
    if [ -s "${relabeled_merged}" ]; then
      # dereplicate to uniques
      vsearch \\
        --sizeout \\
        --threads ${task.cpus} \\
        --derep_fulllength ${relabeled_merged} \\
        --output "${id}_unique.fasta"
    else
      >&2 echo "Merged FASTA is empty. Did your PCR primers match anything?"
      exit 1
    fi
    """
  } else {
    """
    echo "usearch: \$(usearch | head -1 | awk '{print \$2}')" >> settings.yml
    if [ -s "${relabeled_merged}" ]; then
      # dereplicate to uniques
      usearch \\
        -fastx_uniques ${relabeled_merged} \\
        -sizeout \\
        -fastaout "${id}_unique.fasta" \\
        -threads ${task.cpus}
    else
      >&2 echo "Merged FASTA is empty. Did your PCR primers match anything?"
      exit 1
    fi
    """
  }
}

// remove chimera sequences
process remove_chimeras {
  label 'denoiser'
  label 'process_full'

  publishDir "${params.outDir}/zotus", mode: params.publishMode

  input:
    tuple val(id), path(uniques), path(chimera_reference)

  output:
    tuple val(id), path("${id}_chimeras_removed.fasta"), emit: result
    path 'settings.yml'
    path 'chimera_map.tsv'
    path '*_chimera_sequences.fasta'
    path '*_chimeras_denovo.fasta', optional: true
    path '*_chimeras_reference.fasta', optional: true


  script:
  if (params.denoiser == "vsearch") {
    """
    echo "vsearch: \$(vsearch --version 2>&1| head -1 | awk  '{print \$2}' | sed 's/,\$//')" >> settings.yml
    # remove chimeras
    if [ -f "${chimera_reference}" ]; then
      # if we have a valid reference file
      # do reference-based chimera removal in addition to denovo
      vsearch \\
        --threads ${task.cpus} \\
        --uchime3_denovo "${uniques}" \\
        --chimeras "${id}_chimeras_denovo.fasta" \\
        --nonchimeras - |\\
        vsearch \\
          --threads ${task.cpus} \\
          --uchime_ref - \\
          --db "${chimera_reference}" \\
          --nonchimeras "${id}_chimeras_removed.fasta" \\
          --chimeras "${id}_chimeras_reference.fasta" 
    else
      # otherwise just do it denovo
      vsearch \\
        --threads ${task.cpus} \\
        --uchime3_denovo "${uniques}" \\
        --nonchimeras "${id}_chimeras_removed.fasta" \\
        --uchimeout chimera_map.tsv \\
        --chimeras "${id}_chimera_sequences.fasta" 
    fi
    """
  } else {
    """
    echo "usearch: \$(usearch | head -1 | awk '{print \$2}')" >> settings.yml
    # remove chimeras
    usearch -uchime3_denovo "${uniques}" \\
      -uchimeout chimera_map.tsv \\
      -chimeras "${id}_chimera_sequences.fasta" \\
      -nonchimeras "${id}_chimeras_removed.fasta"
    """
  }
}

// denoise to zotus
process denoise {
  label 'denoiser'
  label 'process_full'

  publishDir "${params.outDir}/zotus", mode: params.publishMode

  input:
    tuple val(id), path(sequences)

  output:
    tuple val(id), path("${id}_zotus.fasta"), emit: result
    path 'settings.yml'


  script:
  if (params.denoiser == "vsearch") {
    """
    echo "vsearch: \$(vsearch --version 2>&1| head -1 | awk  '{print \$2}' | sed 's/,\$//')" >> settings.yml
    echo 'min-abundance: ${params.minAbundance}' >> settings.yml
    echo 'alpha: ${params.alpha}' >> settings.yml

    # denoise to zotus
    vsearch \\
      --threads ${task.cpus} \\
      --cluster_unoise "${sequences}" \\
      --centroids "${id}_zotus.fasta" \\
      --minsize ${params.minAbundance}  \\
      --unoise_alpha ${params.alpha} \\
      --relabel Zotu
    """
  } else {
    """
    echo "usearch: \$(usearch | head -1 | awk '{print \$2}')" >> settings.yml
    echo 'min-abundance: ${params.minAbundance}' >> settings.yml
    echo 'alpha: ${params.alpha}' >> settings.yml

    # denoise to zotus
    usearch -unoise3 "${sequences}"  \\
      -zotus "${id}_zotus.fasta" \\
      -threads ${task.cpus} \\
      -minsize ${params.minAbundance} \\
      -unoise_alpha ${params.alpha}
    """
  }
}

// generate zotu table
process generate_sequence_table {
  label 'denoiser'
  label 'process_full'

  publishDir "${params.outDir}/zotus", mode: params.publishMode

  input:
    tuple val(id), path(zotus), path(relabeled_merged)

  output:
    tuple val(id), path('zotu_table.tsv'), emit: result
    path 'zotu_map.tsv'
    path 'settings.yml'


  script:
  if (params.denoiser == "vsearch") {
    """
    echo "vsearch: \$(vsearch --version 2>&1| head -1 | awk  '{print \$2}' | sed 's/,\$//')" >> settings.yml
    echo 'zotu-identity: ${params.zotuIdentity}' >> settings.yml

    # generate zotu table
    vsearch \\
      --threads ${task.cpus} \\
      --usearch_global "${relabeled_merged}" \\
      --db "${zotus}" \\
      --id ${params.zotuIdentity} \\
      --otutabout zotu_table.tsv \\
      --userout zotu_map.tsv \\
      --userfields "query+target" \\
      --top_hits_only
    """
  } else {
    """
    echo "usearch: \$(usearch | head -1 | awk '{print \$2}')" >> settings.yml
    echo 'zotu-identity: ${params.zotuIdentity}' >> settings.yml

    # generate zotu table
    usearch -otutab "${relabeled_merged}" \\
      -id ${params.zotuIdentity} \\
      -threads ${task.cpus} \\
      -zotus "${zotus}" \\
      -otutabout zotu_table.tsv \\
      -mapout zotu_map.tsv
    """
  }
}


// run blast query
process blast {
  label 'blast'
  label 'process_full'

  publishDir { "${params.outDir}/blast" }, mode: params.publishMode, pattern: 'settings.yml'
  publishDir { "${params.outDir}/blast/${db_name}" }, mode: params.publishMode
  publishDir {
    def pid = String.format("%d",(Integer)num(params.percentIdentity ))
    def evalue = String.format("%.3f",num(params.evalue))
    def qcov = String.format("%d",(Integer)num(params.qcov))
    "${params.outDir}/blast/pid${pid}_eval${evalue}_qcov${qcov}_max${params.maxQueryResults}/${db_name}"
  }, mode: params.publishMode

  input:
    tuple path(zotus_fasta), val(db_name), path(db_files), path(taxdb), val(taxids), val(method)

  output:
    tuple val(db_name), path("blast_result.tsv"), emit: result
    path 'settings.yml'

  script:

  // format settings values
  def pid = String.format("%d",(Integer)num(params.percentIdentity))
  def evalue = String.format("%.3f",num(params.evalue))
  def qcov = String.format("%d",(Integer)num(params.qcov))

  // setup and populate the "basic" blast options
  def blast_options = [:]
  // blast_options['task'] = "blastn"
  blast_options['perc_identity'] = params.percentIdentity
  blast_options['evalue'] = params.evalue
  blast_options['qcov_hsp_perc'] = params.qcov
  blast_options['max_target_seqs'] = params.maxQueryResults
  blast_options['task'] = params.blastTask

  // get any --blastn-xxx arguments that may exist
  def blastn_map = task.ext.blastn_map
  // construct -taxids or -negative_taxids argument
  if (taxids && method) {
    taxids = (([]+taxids) - "").join(",")
    def tt = blastn_map[method] ?: ""
    blastn_map[method] = ([taxids,tt] - "").join(",")
  } 

  // collapse them into a single string
  def blast_opt_str = (blast_options + blastn_map)
    .collect { k, v -> v == true ? "-${k}" : "-${k} ${v}" }
    .join(" ")
  """
  # record blast settings
  echo "blastn: \$(blastn -version | head -1 | awk '{print \$NF}')" >> settings.yml
  echo "percent-identity: ${params.percentIdentity}" >> settings.yml
  echo "evalue: ${params.evalue}" >> settings.yml
  echo "qcov: ${params.qcov}" >> settings.yml
  echo "max-query-results: ${params.maxQueryResults}" >> settings.yml
  if [ ${blastn_map.size()} -gt 0 ]; then
    echo "blastn-options:" >> settings.yml
    echo -e "${task.ext.blastn_map.collect { k, v -> "  ${k}: ${v}"}.join("\\n")}" >> settings.yml
  fi

  # set BLASTDB to local working directory
  export BLASTDB=.

  # blast our zotus
  blastn \\
    -db "${db_name}" \\
    -outfmt "6 qseqid sseqid staxid ssciname scomname sskingdom pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" \\
    ${blast_opt_str} \\
    -query ${zotus_fasta} -num_threads ${task.cpus} \\
    > blast_result.tsv
  """
}

// merge split blast results
// it's a process instead of collectFile because it has more than one output directory
process merge_split_blasts {
  label 'shell'
  label 'process_single'

  publishDir { "${params.outDir}/blast/${db_name}" }, mode: params.publishMode
  publishDir {
    def pid = String.format("%d",(Integer)num(params.percentIdentity ))
    def evalue = String.format("%.3f",num(params.evalue))
    def qcov = String.format("%d",(Integer)num(params.qcov))
    "${params.outDir}/blast/pid${pid}_eval${evalue}_qcov${qcov}_max${params.maxQueryResults}/${db_name}"
  }, mode: params.publishMode

  input:
    tuple val(db_name), path('result*.tsv')
  output:
    path('blast_result.tsv')

  script:
  """
  cat result*.tsv > blast_result.tsv
  """
}

// lookup taxids from taxa names
process lookup_blast_taxids {
  label 'shell'
  label 'process_single'

  input:
    tuple val(taxa), path(ncbi_dumps)
  output:
    env(taxids)

  script:
  def begin = "BEGIN { " + ([] + taxa).collect { "spp[\"${it.toLowerCase()}\"] = 1;" }.join(" ") + " }"
  """
  taxids=\$(awk -F '\\t' '${begin} (tolower(\$3) in spp && \$7 == "scientific name") {print \$1}' names.dmp | sort -n | paste -sd,)
  """
}

// merge blast results
// it's a process instead of collectFile because it has more than one output directory
process merge_blast {
  label 'shell'
  label 'process_single'

  publishDir {
    def pid = String.format("%d",(Integer)num(params.percentIdentity ))
    def evalue = String.format("%.3f",num(params.evalue))
    def qcov = String.format("%d",(Integer)num(params.qcov))
    "${params.outDir}/blast/pid${pid}_eval${evalue}_qcov${qcov}_max${params.maxQueryResults}"
  }, mode: params.publishMode
  publishDir { "${params.outDir}/blast" }, mode: params.publishMode

  input:
    path 'staged/??????????.tsv'

  output:
    path 'blast_result_merged.tsv'

  script:
  """
  cat staged/*.tsv | sort -k1,1V > blast_result_merged.tsv
  """
}

// make custom blast database for LULU curation
process lulu_blast {
  label 'blast'
  label 'process_medium'

  input:
    tuple val(key), path(zotus_fasta), path(zotu_table)

  output:
    tuple val(key), path('match_list.txt'), path(zotu_table)

  script:
  """
  # blast zotus against themselves to create the match list LULU needs
  makeblastdb -in ${zotus_fasta} -parse_seqids -dbtype nucl -out ${key}_zotus
  blastn -db ${key}_zotus \\
    -outfmt "6 qseqid sseqid pident" \\
    -out match_list.txt -qcov_hsp_perc 80 \\
    -perc_identity 84 -query ${zotus_fasta} \\
    -num_threads ${task.cpus}
  """
}

// LULU curation
process lulu {
  label 'r'
  label 'process_single'
  label 'process_more_memory'

  publishDir "${params.outDir}/lulu", mode: params.publishMode

  input:
    tuple val(key), path(match_list), path(zotu_table)

  output:
    tuple path("lulu_zotu_table.tsv"), path("lulu_zotu_map.tsv"), path("lulu_result_object.rds"), emit: result
    path 'settings.yml'

  script:
  """
  echo "lulu: \$(Rscript -e 'cat(as.character(packageVersion(\"lulu\")),\"\\n\")')" >> settings.yml
  echo "lulu-min-ratio: ${params.luluMinRatio}" >> settings.yml
  echo "lulu-min-ratio-type: ${params.luluMinRatioType}" >> settings.yml
  echo "lulu-min-match: ${params.luluMinMatch}" >> settings.yml
  echo "lulu-min-rc: ${params.luluMinRc}" >> settings.yml


  lulu.R \\
    -m ${params.luluMinRatio} \\
    -t ${params.luluMinRatioType} \\
    -a ${params.luluMinMatch} \\
    -r ${params.luluMinRc} \\
    ${zotu_table} ${match_list} "lulu_zotu_table.tsv" "lulu_zotu_map.tsv" "lulu_result_object.rds"
  """
}

// assign/collapse taxonomy using the R port of the original python script
process collapse_taxonomy {
  label 'r'
  label 'process_single'
  label 'process_more_memory'

  publishDir {
    "${params.outDir}/taxonomy/lca/qcov${params.lcaQcov}_pid${params.lcaPid}_diff${params.lcaDiff}"
  }, mode: params.publishMode
  publishDir { "${params.outDir}/taxonomy/lca" }

  input:
    tuple path(blast_result), path(dmp), path(lineage)

  output:
    path("lca_taxonomy.tsv"), emit: taxonomy
    path("lca_intermediate.tsv")
    path 'settings.yml'


  script:
  def pf = []
  params.lcaFilterMaxQcov && pf << "--filter-max-qcov"
  params.lcaCaseInsensitive && pf << "--case-insensitive"
  """
  # save settings
  echo "lca-qcov: ${params.lcaQcov}" > settings.yml
  echo "lca-pid: ${params.lcaPid}" >> settings.yml
  echo "lca-diff: ${params.lcaDiff}" >> settings.yml
  echo "lca-filter-max-qcov: ${params.lcaFilterMaxQcov ? 'yes' : 'no'}" >> settings.yml
  echo "lca-taxon-filter: ${params.lcaTaxonFilter}" >> settings.yml
  echo "lca-case-insensitive: ${!params.lcaCaseInsensitive ? 'yes' : 'no'}" >> settings.yml

  collapse_taxonomy.R \\
    --qcov ${params.lcaQcov} \\
    --pid ${params.lcaPid} \\
    --diff ${params.lcaDiff} \\
    --evalue ${params.lcaEvalue} \\
    --merged merged.dmp \\
    --nodes nodes.dmp \\
    --taxid-lineage taxidlineage.dmp \\
    --output lca_taxonomy.tsv \\
    --dropped "${params.dropped}" \\
    --intermediate "lca_intermediate.tsv" \\
    ${pf.join(" ")} \\
    --taxon-filter "${params.lcaTaxonFilter}" \\
    ${blast_result} ${lineage}
  """
}

// run insect classifier model
process insect {
  label 'r'
  label 'process_full'

  publishDir {
    def offs = String.format("%d",(Integer)num(params.insectOffset))
    def thresh = String.format("%.2f",num(params.insectThreshold))
    def minc = String.format("%d",(Integer)num(params.insectMinCount))
    def ping = String.format("%.2f",num(params.insectPing))
    "${params.outDir}/taxonomy/insect/thresh${thresh}_offset${offs}_mincount${minc}_ping${ping}"
  }, mode: params.publishMode
  publishDir { "${params.outDir}/taxonomy/insect" }

  input:
    tuple path(classifier), path(zotus), path(dmp)

  output:
    path('insect_taxonomy.tsv'), emit: taxonomy
    path('insect_model.rds')
    path('settings.yml')

  script:
  def offs = String.format("%d",(Integer)num(params.insectOffset))
  def thresh = String.format("%.2f",num(params.insectThreshold))
  def minc = String.format("%d",(Integer)num(params.insectMinCount))
  def ping = String.format("%.2f",num(params.insectPing))

  """
  # record insect settings
  echo "insect: \$(Rscript -e 'cat(as.character(packageVersion(\"insect\")),\"\\n\")')" >> settings.yml
  echo "insect-offset: ${params.insectOffset}" >> settings.yml
  echo "insect-threshold: ${params.insectThreshold}" >> settings.yml
  echo "insect-min-count: ${params.insectMinCount}" >> settings.yml
  echo "insect-ping: ${params.insectPing}" >> settings.yml

  if [ "${classifier}" != "insect_model.rds" ]; then
    mv ${classifier} insect_model.rds
  fi
  insect.R \\
     --cores ${task.cpus} \\
     --threshold ${params.insectThreshold} \\
     --offset ${params.insectOffset} \\
     --min-count ${params.insectMinCount} \\
     --ping ${params.insectPing} \\
     --lineage rankedlineage.dmp \\
     --output insect_taxonomy.tsv \\
     --merged merged.dmp \\
     ${zotus} insect_model.rds
  """
}

// merge split insect results
// it's a process instead of collectFile because it has more than one output directory
process merge_split_insect {
  label 'shell'
  label 'process_single'

  publishDir {
    def offs = String.format("%d",(Integer)num(params.insectOffset))
    def thresh = String.format("%.2f",num(params.insectThreshold))
    def minc = String.format("%d",(Integer)num(params.insectMinCount))
    def ping = String.format("%.2f",num(params.insectPing))
    "${params.outDir}/taxonomy/insect/thresh${thresh}_offset${offs}_mincount${minc}_ping${ping}"
  }, mode: params.publishMode
  publishDir { "${params.outDir}/taxonomy/insect" }

  input:
    path('insect*.tsv')
  output:
    path('insect_taxonomy.tsv')

  script:
  """
  # preserve header
  head -1 insect1.tsv > insect_taxonomy.tsv
  tail -qn+2 insect*.tsv | sort -k1,1V >> insect_taxonomy.tsv
  """
}

// produce a phyloseq object from pipeline output
process phyloseq {
  label 'r'
  label 'process_single'
  label 'process_more_memory'

  publishDir "${params.outDir}/phyloseq", mode: params.publishMode

  input:
    path(zotu_table)
    path(taxonomy)
    path(metadata)
    path(sequences)

  output:
    path("phyloseq.rds")

  script:
  def opt = []
  params.tree && opt << "--tree"
  params.tree && opt << "--sequences \"${sequences}\""
  params.optimizeTree && opt << "--optimize"
  """
  phyloseq.R \\
    --out phyloseq.rds \\
    ${opt.join(" ")} \\
    "${zotu_table}" "${taxonomy}" "${metadata}"
  """
}

process finalize {
  label 'r'
  label 'process_single'
  label 'process_more_memory'

  publishDir {
    def td = params.standaloneTaxonomy ? 'final/standalone' : 'final'
    "${params.outDir}/${td}"
  }, mode: 'copy'

  input:
    tuple path(zotu_table), path(curated_zotu_table), path(lca_taxonomy), path(insect_taxonomy)

  output:
    path("zotu_table_raw.tsv")
    path("taxonomy.tsv"), emit: taxonomy
    path("zotu_table_final*.tsv")
    path("zotu_table_lca.tsv"), optional: true
    path("zotu_table_insect.tsv"), optional: true

  script:
  def opt = []
  params.abundanceFilter && opt << "--abundance-filter"
  params.rarefy && opt << "--rarefy"
  params.filterMinimum && opt << "--filter-min"
  params.lcaTable && opt << "--lca-table"
  params.insectTable && opt << "--insect-table"

  """
  finalize.R \\
    --filter "${params.taxonFilter}" \\
    --remap "${params.taxonRemap}" \\
    --insect "${insect_taxonomy}" \\
    --lca "${lca_taxonomy}" \\
    --dropped "${params.dropped}" \\
    --controls "${params.controls}" \\
    --control-action "${params.controlAction}" \\
    --control-threshold "${params.controlThreshold}" \\
    --decontam-method "${params.decontamMethod}" \\
    --concentration "${params.dnaConcentration}" \\
    --abundance-threshold "${params.abundanceThreshold}" \\
    --rarefaction-method "${params.rarefactionMethod}" \\
    --permutations "${params.permutations}" \\
    --taxon-priority "${params.taxonPriority}" \\
    --curated "${curated_zotu_table}" \\
    ${opt.join(" ")} \\
    ${zotu_table}
  """
}


// we reuse fastqc/multiqc processes at different steps so they're
// included from an external module
include { fastqc as first_fastqc }    from './modules/modules.nf'
include { fastqc as second_fastqc }   from './modules/modules.nf'
include { multiqc as first_multiqc }  from './modules/modules.nf'
include { multiqc as second_multiqc } from './modules/modules.nf'
include { extract_zip as extract_ncbi_taxonomy } from './modules/modules.nf'
include { extract_targz as extract_ncbi_taxdb } from './modules/modules.nf'

workflow {
  // make sure our arguments are all in order
  check_params()

  def directions = []
  // NCBI archives and the files to extract from them
  def ncbi_taxdumps = ['merged.dmp','nodes.dmp','taxidlineage.dmp','rankedlineage.dmp', 'names.dmp']
  def ncbi_taxdb = 'https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz'
  def ncbi_taxdbs = ['taxdb.bti','taxdb.btd','taxonomy4blast.sqlite3']

  // do standalone taxonomy assignment
  if (params.standaloneTaxonomy) {

    if (!params.noTaxdump) {
      // load and extract NCBI taxonomy
      Channel.fromPath(params.ncbiTaxdump,glob:false) |
        combine(Channel.of(ncbi_taxdumps).toList()) |
        extract_ncbi_taxonomy 

      // collate extracted files into a list channel
      ncbi_dumps = extract_ncbi_taxonomy.out.file |
        toList
    } else {
      ncbi_dumps = Channel.of(ncbi_taxdumps.collect { file(it) })
    }

    // do lca
    if (params.lca) {
      // build blast result channel
      blast_result = Channel.fromPath(params.blastFile, checkIfExists: true)

      blast_result | 
        combine(ncbi_dumps) | 
        combine(Channel.fromPath(params.lcaLineage)) |
        collapse_taxonomy

      // pull out lca table
      collapse_taxonomy.out.taxonomy |
        set { lca_taxonomy }
    } else {
      lca_taxonomy = Channel.fromPath('nofile-lca-taxonomy')
    }

    // do insect
    if (params.insect) {
      // load sequences fasta
      sequences = Channel.fromPath(params.insectSequences, checkIfExists: true)

      // load the classifier model
      if (helper.file_exists(params.insect) || helper.is_url(params.insect)) {
        classifier = Channel.fromPath(params.insect, glob:false)
      } else {
        // download the classifier model if it's one of the supported ones
        // previous sanity checks ensure the model is in our helper map
        def m = params.insect.toLowerCase()
        def url = helper.insect_classifiers[m]
        // glob:false is necessary because the urls have question marks in them
        classifier = Channel.fromPath(url, glob:false)
      }


      // if requested, split query sequences into chunks
      if (params.splitSequences) {
        sequences |
          splitFasta(by: params.splitSequencesBy, file: true) |
          set { query_sequences }
      } else {
        query_sequences = sequences
      }

      // run the insect classification
      classifier |
        combine(query_sequences) |
        combine(ncbi_dumps) |
        insect 

      if (params.splitSequences) {
        insect.out.taxonomy |
          toList |
          merge_split_insect |
          set { insect_taxonomy }
      } else {
        insect_taxonomy = insect.out.taxonomy
      }

    } else {
      insect_taxonomy = Channel.fromPath('nofile-insect-taxonomy')
    }

    // do this part if the sequence table exists
    if (helper.file_exists(params.seqTable)) {
      zotu_table = Channel.fromPath(params.seqTable, checkIfExists: true)
      curated_zotu_table = Channel.fromPath("nofile-curated-zotu-table")

      // run it through finalize
      zotu_table |
        combine(curated_zotu_table) |
        combine(lca_taxonomy) |
        combine(insect_taxonomy) |
        finalize
    }
  } else {
    // save the config file in yaml format
    // TODO: saving the config file to a new file seems to break caching for the blast process?
    // but downstream things that rely on blast are still cached. it makes no sense
    if (params.saveConfig) {
      // default to 'options.yml' in the launch directory
      def config_file = launchDir / "options.yml"
      // if it's a string, use that as the location
      if (params.saveConfig instanceof String) {
        config_file = params.saveConfig
      }
      save_config(config_file.toString())
    }

    // if there isn't an already-demultiplexed FASTA file
    // figure out where the sequence reads are, make sure they're
    // in the right order, and remap sample IDs (if requested)
    if (!helper.file_exists(params.sequences)) {
      if (params.single) {
        // if params.reads is a directory, make it a glob
        def reads_files = params.reads
        if (helper.is_dir(reads_files)) {
          reads_files = file(reads_files) / '*.f*q*'
        }
        // here we load whatever was passed as the --reads option
        // if it's a glob, we get a list of files. if it's just one, we get just one
        // if it's a directory, it's made into a glob to find reads in that directory
        // and we use the basename of the file as a sample ID.
        Channel.fromPath(reads_files, checkIfExists: true) |
          map { [ it.baseName, [it] ] } |
          set { reads }
      } else if (params.paired) {
        // if fwd and rev point to files that exists, just load them directly
        if ( helper.file_exists(params.fwd) && helper.file_exists(params.rev) ) {
          directions = [params.fwd,params.rev]
          Channel.of(params.project) |
            combine(Channel.fromPath(params.fwd,checkIfExists: true)) |
            combine(Channel.fromPath(params.rev,checkIfExists: true)) |
            map { a,b,c -> [a,[b,c]] } |
            set { reads }
         } else {
          // figure out how the reads are to be found and find them
          def pattern = ""
          // if --fwd and --rev are both globs
          if ( params.fwd != "" && params.rev != "" && helper.is_list(file(params.fwd)) && helper.is_list(file(params.rev)) ) {
            // get directory part of fwd and rev globs.
            // file() resolves the glob's fully qualified path
            // and .Parent gets the directory part
            // file() will also get all matches to the glob, so we call unique()
            // to collapse directories, hoping there is only one
            def fwd_path = file(params.fwd).Parent.unique()
            def rev_path = file(params.rev).Parent.unique()

            // make sure globs actually matched something
            if (fwd_path.size() == 0) {
              exit(1,"No files matched by --fwd glob.")
            }
            // make sure globs actually matched something
            if (rev_path.size() == 0) {
              exit(1,"No files matched by --rev glob.")
            }

            // make sure globs only mached a single directory
            if (fwd_path.size() > 1 || rev_path.size() > 1) {
              exit(1,"Files matched by --fwd/--rev globs must reside in single directories.")
            }
            // reduce to first element and convert to string
            fwd_path = fwd_path[0].toString()
            rev_path = rev_path[0].toString()

            // make sure the directory part ends in '/'
            if (fwd_path[-1] != '/') fwd_path += '/'
            if (rev_path[-1] != '/') rev_path += '/'

            // CEB: The strategy here is to idenitify identical and non-identical text in the --fwd and --rev globs to generate
            // a single glob that is compatible with Channel.fromFilePairs

            // Extract the glob part from the provided paths
            // new File() is used because it parses but does not resolve globs
            def fwd_file_pattern = new File(params.fwd).Name
            def rev_file_pattern = new File(params.rev).Name

            // get common prefix
            // for some weird reason, nextflow doesn't like having def and the assignment on the same line
            // when assigning the results of a call to a helper class method, so define it first
            def common_path_prefix = ""
            common_path_prefix = helper.common(fwd_path,rev_path)

            // CEB: Adjust common_prefix to remove the last character if it's where they start to differ
            // MH: I'm not sure why this is needed but it doesn't seem to harm anything so I'm leaving it in
            if (common_path_prefix.size() > 0 && fwd_path.charAt(common_path_prefix.size() - 1) != rev_path.charAt(common_path_prefix.size() - 1)) {
              common_path_prefix = common_path_prefix[0..-2]
            }

            // CEB: Only run the following block if the fwd_path and rev_path are different
            def fwd_path_diff = ""
            def rev_path_diff = ""
            def common_path_suffix = ""
            if (fwd_path != rev_path) {
              // Extract the common path suffix by comparing the strings in reverse
              def fwd_path_reversed = fwd_path.reverse()
              def rev_path_reversed = rev_path.reverse()

              common_path_suffix = helper.common(fwd_path_reversed,rev_path_reversed).reverse()

              // Extract the differing middle part of the path
              fwd_path_diff = fwd_path.substring(common_path_prefix.size(), fwd_path.size() - common_path_suffix.size())
              rev_path_diff = rev_path.substring(common_path_prefix.size(), rev_path.size() - common_path_suffix.size())

            }
            // CEB: Extract the common file prefix
            // for some weird reason, nextflow doesn't like having def and the assignment on the same line
            // when assigning the results of a call to a helper class method, so define it first
            def common_file_prefix = ""
            common_file_prefix = helper.common(fwd_file_pattern,rev_file_pattern)

            // CEB: Adjust common_file_prefix to remove the last character if it's where they start to differ
            // MH: (again, not sure why this is necessary?)
            if (common_file_prefix.size() > 0 && fwd_file_pattern.charAt(common_file_prefix.size() - 1) != rev_file_pattern.charAt(common_file_prefix.size() - 1)) {
              common_file_prefix = common_file_prefix[0..-2]
            }

            // CEB: Extract the common file suffix by comparing the strings in reverse
            def fwd_file_pattern_reversed = fwd_file_pattern.reverse()
            def rev_file_pattern_reversed = rev_file_pattern.reverse()
            // for some weird reason, nextflow doesn't like having def and the assignment on the same line
            // when assigning the results of a call to a helper class method, so define it first
            def common_file_suffix = ""
            common_file_suffix = helper.common(fwd_file_pattern_reversed,rev_file_pattern_reversed).reverse()

            // CEB: Extract the differing middle part of the file pattern
            def fwd_file_diff = ""
            def rev_file_diff = ""
            if (common_file_prefix != common_file_suffix) {
              // MH: check for certain edge cases
              if (common_file_suffix.size() >= fwd_file_pattern.size())
                fwd_file_diff == ""
              else
                fwd_file_diff = fwd_file_pattern.substring(common_file_prefix.size(), fwd_file_pattern.size() - common_file_suffix.size())
              if (common_file_suffix.size() >= rev_file_pattern.size())
                rev_file_diff == ""
              else
                rev_file_diff = rev_file_pattern.substring(common_file_prefix.size(), rev_file_pattern.size() - common_file_suffix.size())
            }

            // MH: if these are both blank, we don't need the '{,}' part
            def fd = (fwd_file_diff + rev_file_diff != "") ? "{${fwd_file_diff},${rev_file_diff}}" : ""
            // MH: a boolean shorthand to check if these are the same
            def suf = common_file_prefix == common_file_suffix

            // CEB: Construct the final pattern
            pattern = fwd_path == rev_path
              ? "${common_path_prefix}${common_file_prefix}{${fwd_file_diff},${rev_file_diff}}${common_file_suffix}"
              : suf ? "${common_path_prefix}{${fwd_path_diff},${rev_path_diff}}${common_path_suffix}${common_file_suffix}"
                : "${common_path_prefix}{${fwd_path_diff},${rev_path_diff}}${common_path_suffix}${common_file_prefix}${fd}${common_file_suffix}"

          // CEB Add support for --reads glob.  Glob must follow rules for NextFlow Channel.fromFilePairs.
          //    Basically, the glob should contain [12] or {R1,R2} or etc... based on my testing
          //    The only way to get away from this requirement is to tell fromFilePairs how many files to expect, or
          //    to write a script that generates a compatible glob from the files returned by --reads glob
          } else if (params.reads != "" && helper.is_list(file(params.reads))) {
            if (file(params.reads).size() > 0)
              pattern = "${params.reads}"
            else exit(1,"No files matched by --reads glob")

          //CEB dirs are specified by --reads, --fwd, --rev, original functionality
          } else if (params.fwd != "" && params.rev != "" && params.reads != "" && helper.is_dir(params.reads + '/' + params.fwd) && helper.is_dir(params.reads + "/" + params.rev)) {
            pattern = "${params.reads}/{${params.fwd},${params.rev}}/*{${params.r1},${params.r2}}*.f*q*"

          //CEB user provides dirs for --fwd and --rev but not --reads, new functionality (borrow code from --fwd --ref globs above)
          } else if (helper.is_dir(params.fwd) && helper.is_dir(params.rev)) {
            // MH: there's a weird business where if the things inside of the {} end with '/',
            // the glob is not matched (even though this works in bash).
            // so we'll strip off any trailing slash
            def f = params.fwd.replaceAll(/\/$/,'')
            def r = params.rev.replaceAll(/\/$/,'')
            pattern = "{${f},${r}}/*{${params.r1},${params.r2}}*.f*q*"

          //CEB dir is specified by --reads; original functionality
          } else if (helper.is_dir(params.reads)) {
            pattern = "${params.reads}/*{${params.r1},${params.r2}}*.f*q*"
          } else {
            exit(1,"Arguments passed to --reads, --fwd, and/or --rev point to directories and/or files that do not exist")
          }

          // Construct reads channel
          // Replace hyphens with underscores in the sample name, because some tools
          // (notably vsearch) will cut on hyphens as a delimiter and potentially cause havok as a result.

          Channel.fromFilePairs(pattern, checkIfExists: true) |
            // make sure we have a key value (project ID)
            map { key,reads -> [ key ?: params.project, reads ] } |
            ifEmpty {
              // bail if we didn't find anything
              exit(1,"No paired reads matched by pattern '${pattern}'. Check command-line options.")
            } |
            set { reads }
        }
      } else {
        println(colors.red("Somehow neither ") + colors.bred("--single") + colors.red(" nor ") + colors.bred("--paired") + colors.red(" were passed and we got to this point"))
        println(colors.red("That should not have happened"))
        exit(1)
      }

      // Attempt to enforce proper read order of read pairs, since
      // Channnel.fromFilePairs will load them alphabetically
      // and they might show up in the wrong order
      if (params.r1 != "" && params.r2 != "") {
        reads |
          map { id, reads -> [ 
            id,
            reads.sort{ a,b ->
              a.baseName =~ /${params.r1}/ && b.baseName =~ /${params.r2}/ ? -1 :
                a.baseName =~ /${params.r2}/ && b.baseName =~ /${params.r1}/ ? 1 : 
                  exit(1,"Unable to determine correct order of sequence read files.")
            } 
          ] } |
          set { reads }
      }

      // remap sample IDs if a sample map was provided
      if (helper.file_exists(params.sampleMap)) {
        Channel.fromPath(params.sampleMap) |
          splitCsv(sep: "\t") |
          map{ it[0] =~ /^#/ ? null : [ it[1..-1].collect{ file(it).baseName }.join("-"), it[0] ] } | 
          set { sample_map }
        reads |
          map{ id, pair -> [ pair.collect{ file(it).baseName }.join("-"), pair ] } |
          join( sample_map ) |
          map{ oldid, pair, newid -> [ newid, pair ] } |
          ifEmpty {
            // bail if we didn't find anything
            exit(1,"Sample re-mapping resulted in an empty dataset. Double check contents of map file `${params.sampleMap}`")
          } |
          set { reads }
      }
    }
      
    // run the dada2 pipeline
    if (params.denoiser == "dada2") {

      def trim = params.barcode || (params.fwdPrimer && params.reversePrimer)
      if (params.demultiplexedBy == "index") {
        if (trim) {
          primers = Channel.of([])
          // if there's a barcode file, assume it's in ngsfilter format
          // and pull unique primer pairs out of the fourth and fifth columns
          // if there's more than one unique set, weird stuff might happen
          if (params.barcode) {
            Channel.fromPath(params.barcode) | 
              splitCsv(sep: "\t") |
              map { !(it[0] =~ /^#/ ) ? [it[3],it[4]] : null } |
              unique |
              set { primers }
          } else {
            // otherwise get the primer sequences from the command line
            primers = Channel.of([params.fwdPrimer,params.reversePrimer])
          }
          trim_primers(reads.combine(primers)) |
            set { reads }
        } else {
          trim_length(reads) |
            map { key, reads -> [ key, helper.is_list(reads) ? reads : [reads] ]} |
            set { reads }
        }
      } else {
        // TODO: do a bunch of demultiplexing and so forth
      }

      // do the quality plots (if requested)
      if (params.plotQualities) { 
        dada_plot_quality_profiles(reads)
      }

      if (!params.plotOnly || !params.plotQualities) {
        // flatten the reads since dada2 works with everything all at once
        samples = reads.collect { it[0] }
        fwd = reads.collect { it[1][0] }
        rev = params.paired ? reads.collect { it[1][1] } : Channel.fromPath('-')
        samples |
          toList |
          combine(fwd.toList()) | 
          combine(rev.toList()) |
          set { to_trim }

        filtered = dada_filter_trim(to_trim).result

        // tuple val(samples), path('*_R1_filtered_trimmed.fastq.gz'), path('*_R2_filtered_trimmed.fastq.gz'), path('fwd.rds'), path('rev.rds'), emit: result
        fwd = filtered.map{ [ it[1], it[3], 'fwd' ] }
        rev = filtered.map{ [ it[2], it[4], 'rev' ] }

        dada_learn_errors(params.paired ? fwd.concat(rev) : fwd) |
          set { errors }
        
        if (params.plotErrors) {
          dada_plot_errors(errors)
        }

        errors | 
          dada_infer_samples |
          set { inferred_samples }
        
        if (params.paired) {
          inferred_samples | 
            collect | 
            map { [ it[0] + it[5], [it[1],it[6]], [it[2],it[7]], [it[3],it[8]] ]} |
            set { to_merge }
          to_merge |
            dada_merge_reads | 
            set { denoised }
        } else {
          inferred_samples |
            map { it[3] } |
            set { denoised } 
        }
        dada_make_asv_table(denoised)
        dada_make_asv_table.out.asv | 
          dada_remove_chimeras

        // tuple path(filter_table), path(filter), path(dada), path(merged), path(asv_table)
        if (params.paired) {
          dada_filter_trim.out.filter |
            combine(to_merge) | 
            map { [ it[0], it[2], it[4] ] } |
            combine(dada_merge_reads.out) |
            combine(dada_remove_chimeras.out.asv | map { it[0] }) |
            set { to_track }
        } else {
          dada_filter_trim.out.filter | 
            combine(dada_filter_trim.out.result | map { it[3] }) |
            combine(dada_infer_samples.out | map { it[3] } ) | 
            combine(Channel.fromPath('-')) |
            combine(dada_remove_chimeras.out.asv | map { it[0] }) | 
            set { to_track }
        }
        dada_track_reads(to_track)
        dada_remove_chimeras.out.asv_table | 
          set { seq_table }
        dada_remove_chimeras.out.fasta | 
          set { sequences }
      }
    } else { // run the u/vsearch pipeline
      if (helper.file_exists(params.sequences)) {
        // we've already demultiplexed and relabeled sequences
        // (presumably from an earlier run of the pipeline), so we can jump to here

        // load the fasta file in usearch/vsearch format
        Channel.fromPath(params.sequences, checkIfExists: true) |
          set { to_dereplicate }
      } else {
        // otherwise do all the various processing bits

        // load barcodes or create a barcode file from primers
        // run them through fix_barcodes if we need to, 
        // which replaces I's with N's in the primer sequences
        if (params.barcode) {
          Channel.fromPath(params.barcode) |
            fix_barcodes |
            set { barcodes }
        } else if (params.fwdPrimer && params.reversePrimer) {
          Channel.of( [params.fwdPrimer.replaceAll(/[Ii]/,'N'), params.reversePrimer.replaceAll(/[Ii]/,'N')]  ) | 
            collectFile { ['barcode.tsv', "marker\tsample\t:\t${it[0]}\t${it[1]}\tseq\n"] } | 
            set { barcodes }
        } else {
          barcodes = Channel.fromPath('-')
        }

        // if the sequences are already demultiplexed by indices, we'll
        // process them separately, including optionally attempting to remove ambiguous indices
        // and ultimately smash them together for vsearch/usearch to do the dereplication
        if (params.demultiplexedBy == "index") {
          // do fastqc/multqc before filtering & merging
          if (params.fastqc) {
            Channel.of("initial") |
              combine(reads) |
              first_fastqc |
              collect(flat: true) |
              toList |
              combine(Channel.of("initial")) |
              first_multiqc
          }

          // run the first part of the pipeline for sequences that have already
          // been demultiplexed by the sequencer
          reads |
            filter_merge 
          filter_merge.out.result |
            set { reads_filtered_merged }

          // do fastqc/multiqc for filtered/merged
          if (params.fastqc) {
            Channel.of("filtered") |
              combine(reads_filtered_merged) |
              second_fastqc |
              collect(flat: true) |
              toList |
              combine(Channel.of("filtered")) |
              second_multiqc
          }

          // remove ambiguous indices, if specified
          if (params.removeAmbiguousIndices) {
            reads_filtered_merged |
              filter_ambiguous_indices |
              set { reads_filtered_merged }
          }

          // with or without the primer mismatch check, do the
          // length filtering and smash results together into one file
          reads_filtered_merged |
            combine(barcodes) |
            set { rfm_barcodes }

          // only run ngsfilter if we have primers
          def trim = params.barcode || (params.fwdPrimer && params.reversePrimer)
          if(trim) {
            rfm_barcodes |
              ngsfilter 
            ngsfilter.out.result |
              set { rfm_barcodes }
          }

          // continue length filtering and whatnot
          rfm_barcodes |
            filter_length 
          filter_length.out.result |
            map { [it[0], it[1]]} |
            // collectFile concatenates multiple possible barcode/primer matches
            collectFile { id, file -> [ "${id}.fastq", file ] } |
            map { [ it.baseName, it ] } |
            // relabel to fasta
            relabel |
            set { relabeled }

          relabeled.result |
            toList | merge_relabeled |
            set { to_dereplicate }

        } else { // demultiplexed by barcode/pool
          // here, reads are demultiplexed by barcodes, so they're either
          // all in one or two fastq files (depending on single vs paired end)
          // or they're pooled such that barcode pairs are reused across index pairs

          // split the input fastqs to increase parallelism, if requested
          if (params.split) {
            if (params.paired) {
              reads |
                // flatten the reads tuple
                map { key, reads -> [key] + reads } |
                // split fastq files
                splitFastq(by: params.splitBy, file: true, pe: true) |
                // rearrange reads tuple so it looks like [key, [R1,R2]]
                map { key, read1, read2 -> [key, [read1,read2]] } |
                set { reads }
            } else {
              // in single-end mode we can just split directly
              reads |
                // flatten the reads tuple
                map { key, reads -> [key] + reads } |
                splitFastq(by: params.splitBy as Integer, file: true) |
                map { key, readfile -> [readfile.baseName, readfile] } |
                set { reads }
            }
          }

          // do initial fastqc step
          if (params.fastqc) {
            Channel.of("initial") |
              combine(reads) |
              first_fastqc
            // if input files are split we'll run them through multiqc
            if (params.split || params.demultiplexedBy == "pool") {
              first_fastqc.out |
                collect(flat: true) |
                toList |
                combine(Channel.of("initial")) |
                first_multiqc
            }
          }

          // do quality filtering and/or paired-end merge
          reads |
            filter_merge 
          filter_merge.out.result |
            set { reads_filtered_merged }

          // post-filtering fastqc step
          if (params.fastqc) {
            Channel.of("filtered") |
              combine(reads_filtered_merged) |
              second_fastqc
            // again run multiqc if split
            if (params.split || params.demultiplexedBy == "pool") {
              second_fastqc.out |
                collect(flat: true) |
                toList |
                combine(Channel.of("filtered")) |
                second_multiqc
            }
          }

          // process pooled barcodes
          if (params.demultiplexedBy == "pool") {
            barcodes |
              // split barcode file into multiples by the first column (key value)
              split_barcodes | flatten |
              // and make it a list of [key, split barcode piece]
              map { [it.baseName.split(/---/)[0], it] } |
              set { barcodes }

            // combines pooled reads with barcode files
            reads_filtered_merged |
              // this gives us a huge mess of combinations and many of them are wrong
              combine(barcodes) |
              // so filter them down to the the ones where the key matches
              filter { key1, f1, key2, f2 -> key1 == key2 } |
              // and make sure they're in a format we expect
              map { key1, f1, key2, f2 -> [key1, f1, f2] } |
              set { reads_barcodes }
          } else {
            // combine reads with barcode file(s)
            reads_filtered_merged |
              combine(barcodes) |
              set { reads_barcodes }
          }

          // run the rest of the pipeline, including demultiplexing, length filtering,
          // splitting, and recombination for dereplication
          reads_barcodes |
            ngsfilter 
          ngsfilter.out.result |
            filter_length 
          filter_length.out.result |
            split_samples |
            // we have to flatten here because we can get results that look like
            // [[sample1,sample2,sample3],[sample1,sample2,sample3]]
            flatten | 
            // collect different files with the same name into concatenated samples
            collectFile |
            // extract sample IDs
            map { [it.baseName, it] } |
            // relabel to fasta
            relabel |
            set { relabeled }

          // collect to single relabeled fasta
          relabeled.result |
            toList | merge_relabeled |
            set { to_dereplicate }
        }
      }
      
      if (!params.preprocessOnly) {
        // dereplicate
        Channel.of(params.project) |
          combine(to_dereplicate) |
          dereplicate |
          set { dereplicated }

        // remove chimeras
        dereplicated.result |
          combine(Channel.fromPath(params.chimeraRef)) |
          remove_chimeras |
          set { chimeras_removed }

        // denoise to zotus
        chimeras_removed.result |
          denoise | 
          set { denoised } 

        // generate zotu table
        denoised.result | 
          combine(to_dereplicate) |
          generate_sequence_table | 
          set { table }

        // set channels used further downstream
        denoised.result | 
          map { it[1] } |
          set { sequences }
        
        table.result | 
          map { it[1] } |
          set { seq_table }

      }
    }
      
    if (!params.preprocessOnly && (!params.plotOnly || !params.plotQualities)) {
      if (params.blastTaxa || params.blastExcludeTaxa || params.insect || params.lca) {
        if (!params.noTaxdump) {
          // load and extract NCBI taxonomy
          Channel.fromPath(params.ncbiTaxdump,glob:false) |
            combine(Channel.of(ncbi_taxdumps).toList()) |
            extract_ncbi_taxonomy 

          // collate extracted files into a list channel
          ncbi_dumps = extract_ncbi_taxonomy.out.file |
            toList
        } else {
          ncbi_dumps = Channel.of([ ncbi_taxdumps.collect { file(it) } ])
        }
      }

      // if requested, split query sequences into chunks
      if (params.splitSequences) {
        sequences |
          splitFasta(by: params.splitSequencesBy, file: true) |
          set { query_sequences }
      } else {
        query_sequences = sequences
      }

      // run blast query, unless skipped
      if (params.blast) {
        // def only works on its own line
        // possibly related to NF issue #804: https://github.com/nextflow-io/nextflow/issues/804

        // get --blast param value(s) as a list with only unique elements
        def blasts = ([] + params.blast).unique()

        // collect list of blast database files, grouped by database name
        Channel.fromPath(blasts) | 
          map { [ it.Name, file("${it}.*") ] } | 
          set { blastdb } 

        // get taxdb if specified on command line
        if (params.blastTaxdb) {
          // stage/download file and extract
          // glob:false required for URLs to work properly
          Channel.fromPath(params.blastTaxdb === true ? ncbi_taxdb : params.blastTaxdb,glob:false) | 
            combine(Channel.of(ncbi_taxdbs).toList()) |
            extract_ncbi_taxdb
          // flatten to list
          extract_ncbi_taxdb.out.file |
            toList |
            set { tdb }
          // combine with blast db channel
          blastdb = blastdb.combine(tdb)
        } else {
          // otherwise just assume taxdb files live under each blast db
          Channel.fromPath(blasts, checkIfExists: false) |
            map { b -> [b.Name, ncbi_taxdbs.collect{ file("${b.Parent}/${it}") } ] } |
            set { tdb }
          blastdb = blastdb.join(tdb)
        }

        // create the blast input channel
        query_sequences | 
          combine(blastdb) | 
          set { blast_input }

        // default taxid filter values are blank
        def blast_filter_method = ''
        blast_taxids = Channel.of(["",""])

        // get taxids to include/exclude if requested
        if (params.blastTaxa || params.blastExcludeTaxa) {
          def taxa = params.blastTaxa ? params.blastTaxa : params.blastExcludeTaxa
          blast_filter_method = params.blastTaxa ? 'taxids' : 'negative_taxids'
          Channel.of(taxa.split(",")).collect().toList() |
            combine(ncbi_dumps) | 
            lookup_blast_taxids |
            toList |
            combine(Channel.of(blast_filter_method)) |
            set { blast_taxids }
        } 
        // run the blast query
        blast(blast_input.combine(blast_taxids))
        blast_result = blast.out.result

        // merge split blast results by database
        if (params.splitSequences) {
          blast_result |
            groupTuple | 
            merge_split_blasts |
            set { blast_result } 
        } else {
          // otherwise just get results
          blast_result | 
            map { it[1] } |
            set { blast_result }
        }

        // merge blast results from different databases
        blast_result |
          collect |
          merge_blast |
          set { blast_result }
      }

      // make lulu blast database and do lulu curation
      if (params.lulu) {
        Channel.of(params.project) | 
          combine(sequences) | 
          combine(seq_table) |
          lulu_blast |
          lulu
      }

      // run the insect classifier, if so desired
      if (params.insect) {
        // load the classifier model
        if (helper.file_exists(params.insect) || helper.is_url(params.insect)) {
          classifier = Channel.fromPath(params.insect, glob:false)
        } else {
          // download the classifier model if it's one of the supported ones
          // previous sanity checks ensure the model is in our helper map
          def m = params.insect.toLowerCase()
          def url = helper.insect_classifiers[m]
          // glob:false is necessary because the urls have question marks in them
          classifier = Channel.fromPath(url, glob:false)
        }

        // run the insect classification
        classifier |
          combine(query_sequences) |
          combine(ncbi_dumps) |
          insect
        if (params.splitSequences) {
          insect.out.taxonomy |
            toList |
            merge_split_insect |
            set { insect_taxonomy }
        } else {
          insect_taxonomy = insect.out.taxonomy
        }
      } else {
        insect_taxonomy = Channel.fromPath("nofile-insect-taxonomy")
      }

      // run taxonomy assignment/collapse script if so requested
      if (params.lca && params.blast) {
        // then we smash it together with the blast results
        // and run the LCA process
        blast_result |
          combine(ncbi_dumps) |
          combine(Channel.fromPath(params.lcaLineage)) |
          collapse_taxonomy
        lca_taxonomy = collapse_taxonomy.out.taxonomy
      } else {
        lca_taxonomy = Channel.fromPath("nofile-lca-taxonomy")
      }

      if (params.lulu) {
        lulu.out.result |
          map { it[0] } |
          set { curated_zotu_table }
      } else {
        curated_zotu_table = Channel.fromPath("nofile-curated-zotu-table")
      }

      // combine and finalize outputs
      if (params.lca || params.insect) {
        seq_table |
          combine(curated_zotu_table) |
          combine(lca_taxonomy) |
          combine(insect_taxonomy) |
          finalize
      }

      // create phyloseq output
      if (params.phyloseq) {
        def physeq = true
        switch (params.taxonomy) {
          case "lca":
            if (params.lca) {
              ph_taxonomy = lca_taxonomy
            } else {
              physeq = false
            }
            break
          case "insect":
            if (params.insect) {
              ph_taxonomy = insect_taxonomy
            } else {
              physeq = false
            }
            break
          case "combined":
            if (params.insect || params.lca) {
              ph_taxonomy = finalize.out.taxonomy
            } else {
              physeq = false
            }
            break
          default:
            if (helper.file_exists(params.taxonomy)) {
              ph_taxonomy = Channel.fromPath(params.taxonomy)
            } else {
              physeq = false
            }
            break
        }
        if (physeq) {
          // construct the phyloseq output
          metadata = Channel.fromPath(params.metadata)
          phyloseq(seq_table,ph_taxonomy,metadata,sequences)
        }
      }
    }
  }
}
