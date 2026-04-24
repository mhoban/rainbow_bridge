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

  // give example of what a demultiplexed FASTA file looks like
  if (params.demuxedExample) {
    helper.demuxed_example()
    exit(0)
  }

  if (params.split && params.demultiplexedBy == "index") {
    println(colors.red("Parameters") + colors.bred(" --split ") + colors.red("and") + colors.bred(" index-based demultiplexing ") + colors.red("are mutually exclusive"))
    exit(1)
  }

  // make sure the right version of single,paired,demultiplexed is passed
  if (!helper.file_exists(params.demuxedFasta) && !params.standaloneTaxonomy && params.single == params.paired) {
    if (!params.single) {
      println(colors.red("One of either ") + colors.bred("--single") + colors.red(" or ") + colors.bred("--paired") + colors.red(" MUST be passed"))
    } else {
      println(colors.red("Only one of either ") + colors.bred("--single") + colors.red(" or ") + colors.bred("--paired") + colors.red(" may be passed"))
    }
    exit(1)
  }

  // check phyloseq params
  if (params.phyloseq) {
    if (!helper.file_exists(params.metadata)) {
      println(colors.yellow("The metadata file you passed to use with phyloseq ('${params.metadata}') does not exist"))
      /* exit(1) */
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
      if (!helper.file_exists(params.insectSequences)) {
        exit(1,colors.red("The supplied FASTA file \"${params.insectSequences}\" does not exist"))
      }
    }
  }

  if (params.sampleMap != "" && !helper.file_exists(params.sampleMap)) {
    println(colors.red("The supplied sample map file ${params.sampleMap} does not exist"))
    exit(1)
  }

  // check to make sure denoiser is a valid input
  if (!(params.denoiser in ['usearch','vsearch'])) {
    exit(1,colors.bred("--denoiser") + colors.red(" must be either 'vsearch' or 'usearch'"))
  }

  // sanity check, blast database
  if (params.blast) {

    // if --blast-taxdb is passed, check that it's a .tar.gz archive
    if (params.blastTaxdb && !(params.blastTaxdb =~ /(?i)\.tar\.gz$/)) {
      println(colors.bred("--blast-taxdb") + colors.red(" must be a .tar.gz archive"))
      exit(1)
    }

    // make --blast-db param into a list, if it isn't
    def blasts = params.blastDb
    if (!helper.is_list(blasts))
      blasts = [blasts]

    // get unique vals
    blasts = blasts.unique(false)

    // make sure we've got at least one db
    if (!blasts.size()) {
      println(colors.red("You must pass at least one value to --blast-db"))
      exit(1)
    } else {
      // make sure all dbs exist
      blasts.each {
        if (!file("${it}.ndb").exists()) {
          println(colors.red("Could not find BLAST database '${it}'. Please provide the path to an existing blast database."))
          if (it =~ /~/) {
            println(colors.yellow("The BLAST database '${it}' contains a tilde ('~') that was not expanded by the shell. Try entering an absolute path."))
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
    if (params.lcaLineage && !helper.file_exists(params.lcaLineage)) {
      exit(1,colors.red("The supplied lineage file \"${params.lcaLineage}\" does not exist"))
    }
  }

  // make sure insect parameter is valid: either a file or one of the pretrained models
  if (params.insect) {
    if (!helper.insect_classifiers.containsKey(params.insect.toLowerCase())) {
      if (!helper.file_exists(params.insect)) {
        println(colors.red("Value passed to ") + colors.bred("--insect") + colors.red(" must be one of the supported builtins or an RDS file"))
        println(colors.red("containing a trained insect classifier model."))
        println(colors.red("See rainbow_bridge.nf ") + colors.bred("--help") + colors.red(" for supported builtin models"))
        exit(1)
      }
    }
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
    echo 'single: true' > settings.yml
    echo 'paired: false' >> settings.yml
    echo 'min-quality: ${params.minQuality}' >> settings.yml
    echo 'max-quality: ${params.maxQuality}' >> settings.yml
    echo 'mate-separator: ${params.mateSeparator}' >> settings.yml

    AdapterRemoval --threads ${task.cpus} --file1 ${reads} \\
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
    echo 'single: false' > settings.yml
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
  echo 'primer-mismatch: ${params.primerMismatch}' > settings.yml

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
  echo 'min-len: ${params.minLen}' > settings.yml

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
    echo 'denoiser: vsearch' > settings.yml

    vsearch --threads ${task.cpus} --fastq_qmax ${params.maxQuality} --fastx_filter ${fastq} --relabel "${key}." --label_suffix ";sample=${key}" --fastaout - | \\
      awk '/^>/ {print;} !/^>/ {print(toupper(\$0))}' > "${key}_relabeled.fasta"
    """
  } else {
    """
    echo 'denoiser: usearch' > settings.yml

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

// dereplication, chimera removal, zOTU table generation
process dereplicate {
  label 'denoiser'
  label 'process_full'

  publishDir "${params.outDir}/zotus", mode: params.publishMode

  input:
    tuple val(id), path(relabeled_merged), path(chimera_reference)

  output:
    tuple val(id), path("${id}_unique.fasta"), path("${id}_zotus.fasta"), path("zotu_table.tsv"), emit: result
    path 'settings.yml'
    path 'zotu_map.tsv'
    path 'chimera_map.tsv'
    path '*_chimera_sequences.fasta'
    path '*_chimeras_denovo.fasta', optional: true
    path '*_chimeras_reference.fasta', optional: true

  script:
  if (params.denoiser == "vsearch") {
    """
    echo 'min-abundance: ${params.minAbundance}' > settings.yml
    echo 'alpha: ${params.alpha}' >> settings.yml
    echo 'zotu-identity: ${params.zotuIdentity}' >> settings.yml

    if [ -s "${relabeled_merged}" ]; then
      # dereplicate to uniques
      vsearch \\
        --sizeout \\
        --threads ${task.cpus} \\
        --derep_fulllength ${relabeled_merged} \\
        --output "${id}_unique.fasta"

      # remove chimeras
      if [ -f "${chimera_reference}" ]; then
        # if we have a valid reference file
        # do reference-based chimera removal in addition to denovo
        vsearch \\
          --threads ${task.cpus} \\
          --uchime3_denovo "${id}_unique.fasta" \\
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
          --uchime3_denovo "${id}_unique.fasta" \\
          --nonchimeras "${id}_chimeras_removed.fasta" \\
          --uchimeout chimera_map.tsv \\
          --chimeras "${id}_chimera_sequences.fasta" 
      fi

      # denoise to zotus
      vsearch \\
        --threads ${task.cpus} \\
        --cluster_unoise "${id}_chimeras_removed.fasta" \\
        --centroids "${id}_zotus.fasta" \\
        --minsize ${params.minAbundance}  \\
        --unoise_alpha ${params.alpha} \\
        --relabel Zotu

      # generate zotu table
      vsearch \\
        --threads ${task.cpus} \\
        --usearch_global ${relabeled_merged} \\
        --db "${id}_zotus.fasta" \\
        --id ${params.zotuIdentity} \\
        --otutabout zotu_table.tsv \\
        --userout zotu_map.tsv \\
        --userfields "query+target" \\
        --top_hits_only
    else
      >&2 echo "Merged FASTA is empty. Did your PCR primers match anything?"
      exit 1
    fi
    """
  } else {
    """
    echo 'min-abundance: ${params.minAbundance}' > settings.yml
    echo 'alpha: ${params.alpha}' >> settings.yml
    echo 'zotu-identity: ${params.zotuIdentity}' >> settings.yml

    if [ -s "${relabeled_merged}" ]; then
      # dereplicate to uniques
      usearch \\
        -fastx_uniques ${relabeled_merged} \\
        -sizeout \\
        -fastaout "${id}_unique.fasta" \\
        -threads ${task.cpus}

      # remove chimeras
      usearch -uchime3_denovo "${id}_unique.fasta" \\
        -uchimeout chimera_map.tsv \\
        -chimeras "${id}_chimera_sequences.fasta" \\
        -nonchimeras "${id}_chimeras_removed.fasta"

      # denoise to zotus
      usearch -unoise3 "${id}_unique.fasta"  \\
        -zotus "${id}_zotus.fasta" \\
        -threads ${task.cpus} \\
        -tabbedout zotu_map.tsv \\
        -minsize ${params.minAbundance} \\
        -unoise_alpha ${params.alpha}

      # generate zotu table
      usearch -otutab ${relabeled_merged} \\
        -id ${params.zotuIdentity} \\
        -threads ${task.cpus} \\
        -zotus ${id}_zotus.fasta \\
        -otutabout zotu_table.tsv \\
        -mapout zotu_map.tsv
    else
      >&2 echo "Merged FASTA is empty. Did your PCR primers match anything?"
      exit 1
    fi
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
    tuple path(zotus_fasta), val(db_name), path(db_files), path(taxdb)
    val taxids

  output:
    path("blast_result.tsv"), emit: result
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
  // blast_options['best_hit_score_edge'] = 0.05
  // blast_options['best_hit_overhang'] = 0.25

  // collapse them into a single string
  def blast_opt_str = blast_options
    .collect { k, v -> v == true ? "-${k}" : "-${k} ${v}" }
    .join(" ")

  def blastn_map = task.ext.blastn_map
  if (taxids instanceof Collection) {
    taxids = taxids.findAll { it != "" }
    taxids = taxids.join(",")
  }
  if (taxids) {
    def tt = blastn_map.containsKey('taxids') ? blastn_map['taxids'] : ""
    blastn_map['taxids'] = ([taxids,tt] - "").join(",")
  }

  def blastn_args = blastn_map
    .collect { k, v -> v == true ? "-${k}" : "-${k} ${v}" }
    .join(" ") 
  """
  # record blast settings
  echo "percent-identity: ${params.percentIdentity}" > settings.yml
  echo "evalue: ${params.evalue}" >> settings.yml
  echo "qcov: ${params.qcov}" >> settings.yml
  echo "max-query-results: ${params.maxQueryResults}" >> settings.yml
  if [ -n "${blastn_args}" ]; then
    echo "blastn-options:" >> settings.yml
    echo -e "${task.ext.blastn_map.collect { k, v -> "  ${k}: ${v}"}.join("\\n")}" >> settings.yml
  fi

  # set BLASTDB to local working directory
  export BLASTDB=.

  # blast our zotus
  blastn \\
    -db "${db_name}" \\
    -outfmt "6 qseqid sseqid staxid ssciname scomname sskingdom pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" \\
    ${blast_opt_str} ${blastn_args} \\
    -query ${zotus_fasta} -num_threads ${task.cpus} \\
    > blast_result.tsv
  """
}

// lookup taxids from taxa names
process lookup_blast_taxids {
  // label 'r'
  label 'shell'
  label 'process_single'

  input:
    tuple val(taxa), path(ncbi_dumps)
  output:
    env(taxids)

  script:
  if (!(taxa instanceof Collection)) {
    taxa = [taxa]
  }
  def begin = "BEGIN { " + taxa.collect { "spp[\"${it.toLowerCase()}\"] = 1;" }.join(" ") + " }"
  """
  taxids=\$(awk -F '\\t' '${begin} (tolower(\$3) in spp && \$7 == "scientific name") {print \$1}' names.dmp | sort -n | paste -sd,)
  """
}

// merge blast results
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
    path 'staged/*'

  output:
    path 'blast_result_merged.tsv'

  script:
  """
  cat staged/* > blast_result_merged.tsv
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
  echo "lulu-min-ratio: ${params.luluMinRatio}" > settings.yml
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
    tuple path(blast_result), path(dmp)

  output:
    path("lca_taxonomy.tsv"), emit: taxonomy
    path("lca_intermediate.tsv")
    path 'settings.yml'


  script:
  def pf = []
  params.lcaFilterMaxQcov && pf << "--filter-max-qcov"
  params.lcaCaseInsensitive && pf << "--case-insensitive"
  def lineage = params.lcaLineage ?: 'rankedlineage.dmp'
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
  echo "insect-offset: ${params.insectOffset}" > settings.yml
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
  // files to extract from ncbi archives
  def ncbi_taxdumps = ['merged.dmp','nodes.dmp','taxidlineage.dmp','rankedlineage.dmp', 'names.dmp']
  def ncbi_taxdbs = ['taxdb.bti','taxdb.btd','taxonomy4blast.sqlite3']

  // do standalone taxonomy assignment
  if (params.standaloneTaxonomy) {

    // load and extract NCBI taxonomy
    Channel.fromPath(params.ncbiTaxdump,glob:false) |
      combine(Channel.of(ncbi_taxdumps).toList()) |
      extract_ncbi_taxonomy 

    // collate extracted files into a list channel
    ncbi_dumps = extract_ncbi_taxonomy.out.file |
      toList

    // do lca
    if (params.lca) {
      // build blast result channel
      blast_result = Channel.fromPath(params.blastFile, checkIfExists: true)

      blast_result | 
        combine(ncbi_dumps) | 
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
      zotus = Channel.fromPath(params.insectSequences, checkIfExists: true)

      // load the classifier model
      if (helper.file_exists(params.insect)) {
        classifier = Channel.fromPath(params.insect)
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
        combine(zotus) |
        combine(ncbi_dumps) |
        insect 
      insect_taxonomy = insect.out.taxonomy
    } else {
      insect_taxonomy = Channel.fromPath('nofile-insect-taxonomy')
    }

    // do this part if the zotu table exists
    if (helper.file_exists(params.zotuTable)) {
      zotu_table = Channel.fromPath(params.zotuTable, checkIfExists: true)
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

    if (!helper.file_exists(params.demuxedFasta)) {
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
          map{ [ it[1..-1].collect{ file(it).baseName }.join("-"), it[0] ] } | 
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

      // load barcodes
      // throw an error if the file(s) are bad, but only if we're not skipping the step that needs them.
      // we don't check in check_params because params.barcode could be a wildcard, which is trickier to check cleanly
      // and is handled automatically by fromPath
      // also, we run it through fix_barcodes, which replaces I's with N's in the primer sequences
      Channel.fromPath(params.barcode, checkIfExists: !params.noPcr) |
        fix_barcodes |
        set { barcodes }

      // if the sequences are already demultiplexed by illumina, we'll
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
        if(!params.noPcr) {
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

      } else { // demultiplexed by barcode/combined
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
              splitFastq(by: params.splitBy, file: true) |
              map { key, readfile -> [key, readfile] } |
              set { reads }
          }
        }

        // do initial fastqc step
        if (params.fastqc) {
          Channel.of("initial") |
            combine(reads) |
            first_fastqc
          // if input files are split we'll run them through multiqc
          if (params.split || params.demultiplexedBy == "combined") {
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
          if (params.split || params.demultiplexedBy == "combined") {
            second_fastqc.out |
              collect(flat: true) |
              toList |
              combine(Channel.of("filtered")) |
              second_multiqc
          }
        }

        // process pooled barcodes
        if (params.demultiplexedBy == "combined") {
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
    } else {
      // here we've already demultiplexed and relabeled sequences
      // (presumably from an earlier run of the pipeline), so we can jump to here

      // load the fasta file in usearch/vsearch format
      Channel.fromPath(params.demuxedFasta, checkIfExists: true) |
        set { to_dereplicate }
    }

    // build the input channel, run dereplication, and set to a channel we can use again
    Channel.of(params.project) |
      combine(to_dereplicate) |
      combine(Channel.fromPath(params.chimeraRef)) |
      dereplicate |
      set { dereplicated }

    dereplicated.result |
      set { dereplicated }

    // load and extract NCBI taxonomy
    Channel.fromPath(params.ncbiTaxdump,glob:false) |
      combine(Channel.of(ncbi_taxdumps).toList()) |
      extract_ncbi_taxonomy 
    // collate extracted files into a list channel
    ncbi_dumps = extract_ncbi_taxonomy.out.file |
      toList

    // run blast query, unless skipped
    if (params.blast) {
      // def only works on its own line
      // possibly related to NF issue #804: https://github.com/nextflow-io/nextflow/issues/804

      // make --blast-db value a list, if it's not already
      def blasts = params.blastDb
      if (!helper.is_list(blasts))
        blasts = [blasts]

      // get unique blast dbs
      blasts = blasts.unique(false)

      // wildcard to capture blast database files
      def wildcard = "{.n*,.[0-9]*.n*}"

      // collect list of files within blast databases
      // and group them by blast db names
      Channel.fromPath(blasts) | 
        map { [ it.Name, file("${it}*") ] } | 
        set { blastdb }

      if (!helper.file_exists(params.lcaLineage)) {
        // make channel for taxdb files (whether or not they actually exist)

        // get taxdb if specified on command line
        if (params.blastTaxdb) {
          // stage/download file and extract
          // glob:false required for URLs to work properly
          Channel.fromPath(params.blastTaxdb,glob:false) | 
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
      } else {
        Channel.value( [ ncbi_taxdbs.collect{ file(it) } ] ) |
          set { tdb }
        blastdb = blastdb.combine(tdb)
      }

      // TODO: there's unmatched braces or something somewhere

      // create the blast input channel
      dereplicated |
        map { sid, uniques, zotus, zotutable -> zotus } |
        combine(blastdb) |
        set { blast_input }

      // lookup filter taxids if necessary
      if (params.blastTaxonFilter) {
        Channel.of(params.blastTaxonFilter.split(",")).collect().toList() |
          combine(ncbi_dumps) | 
          lookup_blast_taxids |
          toList |
          set { blast_taxids }
      } else {
        blast_taxids = Channel.of("")
      }

      // run the blast query
      blast(blast_input,blast_taxids)

      // merge blast results from different databases
      blast.out.result |
        collect |
        merge_blast |
        set { blast_result }
    }

    // grab the zotu table from our dereplication step
    dereplicated |
      map { sid, uniques, zotus, zotutable -> zotutable } |
      set { zotu_table }

    // make lulu blast database and do lulu curation
    if (params.lulu) {
      dereplicated |
        // get zotus and sample id
        map { sid, uniques, zotus, zotutable -> [sid,zotus,zotutable] } |
        lulu_blast |
        lulu
    }

    // run the insect classifier, if so desired
    if (params.insect) {
      // dereplicate returns a tuple, but we only need the zotus fasta
      dereplicated |
        map { sid, uniques, zotus, zotutable -> [zotus] } |
        set { zotus }

      // load the classifier model
      if (helper.file_exists(params.insect)) {
        classifier = Channel.fromPath(params.insect)
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
        combine(zotus) |
        combine(ncbi_dumps) |
        insect
      insect_taxonomy = insect.out.taxonomy
    } else {
      insect_taxonomy = Channel.fromPath("nofile-insect-taxonomy")
    }

    // run taxonomy assignment/collapse script if so requested
    if (params.lca && params.blast) {
      // then we smash it together with the blast results
      // and run the LCA process
      blast_result |
        combine(ncbi_dumps) |
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
      zotu_table |
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
        seqs = dereplicated.map { sid, uniques, zotus, zotutable -> zotus }
        phyloseq(zotu_table,ph_taxonomy,metadata,seqs)
      }
    }
  }
}
