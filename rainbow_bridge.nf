#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// pull in the helper class
import helper
import colors

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
        if (!params.collapseTaxonomy) {
          println(colors.yellow("You passed --phyloseq with 'lca' as the taxonomy option, but LCA has not been run."))
          println(colors.yellow("Did you forget the --collapse-taxonomy option?"))
        }
        break
      case 'insect':
        if (!params.insect) {
          println(colors.yellow("You passed --phyloseq with 'insect' as the taxonomy option, but insect has not been run."))
          println(colors.yellow("Did you forget the --insect option?"))
        }
        break
      case 'combined':
        if (!params.collapseTaxonomy && !params.insect) {
          println(colors.yellow("Note: phyloseq generation requires one of --insect, --collapse-taxonomy, or a custom taxonomy table."))
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
    if (!helper.file_exists(params.blastFile)) {
      println(colors.red("The supplied blast result table \"${params.blastFile}\" does not exist"))
      exit(1)
    }
  }

  if (params.sampleMap != "" && !helper.file_exists(params.sampleMap)) {
    println(colors.red("The supplied sample map file ${params.sampleMap} does not exist"))
    exit(1)
  }

  // the --vsearch param is no longer supported, but try to catch if someone still uses it
  if (params.containsKey('vsearch')) {
    println(colors.yellow("FYI: vsearch is now used as the default denoiser and the ") + colors.byellow("--vsearch") + colors.yellow(" option is ignored.\nIf you want to use usearch, pass --denoiser usearch."))
  }

  // if denoiser is an executable, treat it as such
  // otherwise check to make sure it's a valid input
  if (!params.execDenoiser) {
    if (!(params.denoiser in ['usearch','usearch32','vsearch'])) {
      println(colors.bred("--denoiser") + colors.red(" must be either 'usearch', 'usearch32', 'vsearch', or a path to an executable (e.g., /opt/sw/bin/usearch64)"))
      exit(1)
    }
  }

  // sanity check, blast database
  if (params.blast) {

    // if --blast-taxdb is passed, check that it's a directory and that files exist
    if (params.blastTaxdb != "nofile-blast-taxdb") {
      def ok = false
      // check for directory
      if (helper.is_dir(params.blastTaxdb)) {
        // make sure we have two files matching our search pattern
        if( file(params.blastTaxdb).list().findAll{it =~ /taxdb\.bt[id]/}.size() == 2) {
          ok = true
        }
      }
      // otherwise give an error message and bail
      if (!ok) {
        println(colors.bred("--blast-taxdb") + colors.red(" must be a directory containing the NCBI BLAST taxonomy files (taxdb.btd, taxdb.bti)"))
        exit(1)
      }
    }

    // get blast environment variable
    def bdb = helper.get_env("FLOW_BLAST")

    // make --blast-db param into a list, if it isn't
    def blasts = params.blastDb
    if (!helper.is_list(blasts))
      blasts = [blasts]

    // if $FLOW_BLAST was set and we're not ignoring it, add it to the front of the list
    if (bdb != "" && !params.ignoreBlastEnv)
      blasts = blasts.plus(0,bdb)

    // get unique vals
    blasts = blasts.unique(false)

    // make sure we've got at least one db
    if (!blasts.size()) {
      println(colors.red("You must pass at least one value to --blast-db or FLOW_BLAST must point to a valid BLAST database"))
      exit(1)
    } else {
      // make sure all dbs exist
      blasts.each {
        if (!file("${it}.ndb").exists()) {
          println(colors.red("Could not find BLAST database '${it}'. Please provide the path to an existing blast database."))
          if (it =~ /~/) {
            println(colors.yellow("The BLAST database '${it}' contains a '~' that was not expanded by the shell. Try entering an absolute path."))
          }
          exit(1)
        }
      }
    }
  } else {
    if (params.collapseTaxonomy) {
      println(colors.red("--collapse-taxonomy requires the --blast option."))
      exit(1)
    }
  }

  // check LCA options
  if (params.collapseTaxonomy) {
    if (params.lcaDiff <= 0) {
      println(colors.red("--lca-diff argument must be a number greater than zero."))
      exit(1)
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

// parse read directions from a glob
// e.g., file*{R1,R2} -> [R1,R2]
def parse_directions(reads) {
  if (m = reads =~ /\{([^,]+),([^}]+)\}/) {
    return [m.group(1),m.group(2)]
  } else if (m = reads =~ /\[(\w)(\w)\]/) {
    return [m.group(1),m.group(2)]
  } else {
    return []
  }
}


// trim and (where relevant) merge paired-end reads
process filter_merge {
  label 'adapterRemoval'
  label 'process_medium'

  publishDir "${params.preDir}/trim_merge", mode: params.publishMode

  input:
    tuple val(key), path(reads)

  output:
    tuple val(key), path('*_trimmed_merged.fastq')

  script:
  if( reads instanceof Path ) {
    // single end
    """
    AdapterRemoval --threads ${task.cpus} --file1 ${reads} \
      --trimns --trimqualities \
      --minquality ${params.minQuality} \
      --qualitymax ${params.maxQuality} \
      --mate-separator ${params.mateSeparator} \
      --basename ${key}

    mv ${key}.truncated ${key}_trimmed_merged.fastq
    """
  } else {
    // if reads are paired-end then merge
    """
    AdapterRemoval --threads ${task.cpus} --file1 ${reads[0]} --file2 ${reads[1]} \
      --collapse --trimns --trimqualities \
      --minquality $params.minQuality \
      --qualitymax ${params.maxQuality} \
      --minalignmentlength ${params.minAlignLen} \
      --mate-separator ${params.mateSeparator} \
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
    path("${barcodes.BaseName}_fixed.${barcodes.Extension}")

  script:
  """
  if [[ -e "${barcodes}" ]]; then
    fix_barcode.awk "${barcodes}" > "${barcodes.BaseName}_fixed.${barcodes.Extension}"
  else
    touch "${barcodes.BaseName}_fixed.${barcodes.Extension}"
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
    tuple val(key), path("*_annotated.fastq"), val("${barcode.baseName}")

  script:
  """
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
    tuple val(key), path('*_length_filtered.fastq'), val(barcode)

  script:
  """
  obigrep --uppercase -l ${params.minLen} "${fastq}" > "${key}_length_filtered.fastq"
  """
}

// for non-demultiplexed runs, split the annotated reads file by samples
process split_samples {
  label 'obitools'
  label 'process_single'

  publishDir "${params.preDir}/split_samples", mode: params.publishMode

  input:
    tuple val(key), path(fastq), val(barcode)

  output:
    path("__split__*.fastq"), optional: true

  script:
  """
  obisplit --uppercase -p "__split__" -t sample -u "${key}_split_orphans.fastq" ${fastq}
  """
}

// relabel files for dereplication
process relabel {
  label 'denoiser'
  label 'process_low'

  publishDir "${params.preDir}/relabeled", mode: params.publishMode

  input:
    tuple val(key), path(fastq, name: 'input-????.fastq')
  output:
    path('*_relabeled.fasta'), optional: true, emit: result
    path 'settings.txt'


  script:

  if (params.denoiser == "vsearch") {
    def combined = "<(cat input-*.fastq)"
    """
    echo "denoiser: vsearch" > settings.txt
    # this may or may not be necessary anymore, but it seems like a good sanity check
    # since this will fail on empty files
    vsearch --threads ${task.cpus} --fastq_qmax ${params.maxQuality} --fastx_filter ${combined} --relabel "${key}." --fastaout - | \
      awk '/^>/ {print;} !/^>/ {print(toupper(\$0))}' > "${key}_relabeled.fasta"
    """
  } else {
    // combined = fastq.collect{ it.baseName }.join('_') + "_combined.fastq"
    def combined = "combined.fastq"
    def denoiser = params.execDenoiser ? params.denoiser : 'usearch'
    """
    echo "denoiser: ${params.denoiser}" > settings.txt
    cat input-*.fastq > ${combined}
    # we have to convert everything to uppercase because obisplit --uppercase is broken
    ${denoiser} -fastx_relabel ${combined} -prefix "${key}." -fastaout /dev/stdout | \
      tail -n+7 | \
      awk '/^>/ {print;} !/^>/ {print(toupper(\$0))}' > "${key}"_relabeled.fasta
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
  label 'process_high'

  publishDir "${params.outDir}/zotus", mode: params.publishMode

  input:
    tuple val(id), path(relabeled_merged)

  output:
    tuple val(id), path("${id}_unique.fasta"), path("${id}_zotus.fasta"), path("zotu_table.tsv"), emit: result
    path 'settings.txt'
    path 'zotu_map.tsv'

  script:
  if (params.denoiser == "vsearch") {
    """
    echo "denoiser: vsearch" > settings.txt
    echo "minimum sequence abundance: ${params.minAbundance}" >> settings.txt
    echo "alpha: ${params.alpha}" >> settings.txt
    echo "fractional identity: ${params.zotuIdentity}" >> settings.txt
    # steps:
    # 1. get unique sequence variants
    # 2. run denoising algorithm
    # 3. get rid of chimeras
    # 4. match original sequences to zotus by 97% identity
    if [ -s "${relabeled_merged}" ]; then
      vsearch \
        --threads ${task.cpus} --fastq_qmax ${params.maxQuality} \
        --derep_fulllength ${relabeled_merged} --sizeout \
        --output "${id}_unique.fasta"
      vsearch \
        --threads ${task.cpus} --fastq_qmax ${params.maxQuality} \
        --cluster_unoise "${id}_unique.fasta" --centroids "${id}_centroids.fasta" \
        --minsize ${params.minAbundance} --unoise_alpha ${params.alpha}
      vsearch \
        --threads ${task.cpus} --fastq_qmax ${params.maxQuality} \
        --uchime3_denovo "${id}_centroids.fasta" --nonchimeras "${id}_zotus.fasta" \
        --relabel Zotu
      vsearch \
        --threads ${task.cpus} --fastq_qmax ${params.maxQuality} \
        --usearch_global ${relabeled_merged} --db "${id}_zotus.fasta" \
        --id ${params.zotuIdentity} --otutabout zotu_table.tsv \
        --userout zotu_map.tsv --userfields "query+target" \
        --top_hits_only
    else
      >&2 echo "Merged FASTA is empty. Did your PCR primers match anything?"
      exit 1
    fi
    """
  } else {
    def denoiser = params.execDenoiser ? params.denoiser : 'usearch'
    """
    echo "denoiser: ${denoiser}" > settings.txt
    echo "minimum sequence abundance: ${params.minAbundance}" >> settings.txt
    echo "alpha: ${params.alpha}" >> settings.txt
    echo "fractional identity: ${params.zotuIdentity}" >> settings.txt
    # steps:
    # 1. get unique sequences
    # 2. run denoising & chimera removal
    # 3. generate zotu table
    if [ -s "${relabeled_merged}" ]; then
      ${denoiser} -fastx_uniques ${relabeled_merged} \
        -sizeout -fastaout "${id}_unique.fasta"
      ${denoiser} -unoise3 "${id}_unique.fasta"  -zotus "${id}_zotus.fasta" \
        -tabbedout "${id}_unique_unoise3.txt" -minsize ${params.minAbundance} \
        -unoise_alpha ${params.alpha}
      ${denoiser} -otutab ${relabeled_merged} -id ${params.zotuIdentity} \
        -zotus ${id}_zotus.fasta -otutabout zotu_table.tsv -mapout zotu_map.tsv
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
  label 'all_cpus'

  publishDir {
    def pid = String.format("%d",(Integer)num(params.percentIdentity ))
    def evalue = String.format("%.3f",num(params.evalue))
    def qcov = String.format("%d",(Integer)num(params.qcov))
    return "${params.outDir}/blast/pid${pid}_eval${evalue}_qcov${qcov}_max${params.maxQueryResults}/${db_name}"
  }, mode: params.publishMode

  input:
    tuple path(zotus_fasta), val(db_name), path(db_files), path(taxdb)

  output:
    path("blast_result.tsv"), emit: result
    path 'blast_settings.txt'

  script:

  // format settings values
  def pid = String.format("%d",(Integer)num(params.percentIdentity))
  def evalue = String.format("%.3f",num(params.evalue))
  def qcov = String.format("%d",(Integer)num(params.qcov))

  // setup and populate the "basic" blast options
  def blast_options = [:]
  blast_options['task'] = "blastn"
  blast_options['perc_identity'] = params.percentIdentity
  blast_options['evalue'] = params.evalue
  blast_options['qcov_hsp_perc'] = params.qcov
  blast_options['max_target_seqs'] = params.maxQueryResults
  blast_options['best_hit_score_edge'] = 0.05
  blast_options['best_hit_overhang'] = 0.25

  // get extra blast options passed on the command line as --blastn-*
  def extra_options = params
    .findAll { it.key =~ /^blastn-/ }
    .collectEntries { k, v -> [k.tokenize('-')[1],v] }

  // merge blast options with any extra options
  blast_options = blast_options << extra_options
  // collapse them into a single string
  def blast_opt_str = blast_options
    .collect { k, v -> v == true ? "-${k}" : "-${k} ${v}" }
    .join(" ")
  """
  # record blast settings
  echo "Percent identity: ${pid}" > blast_settings.txt
  echo "e-value: ${evalue}" >> blast_settings.txt
  echo "Query qoverage: ${qcov}" >> blast_settings.txt
  echo "Max. target sequences: ${params.maxQueryResults}" >> blast_settings.txt

  # set BLASTDB to local working directory
  export BLASTDB=.

  # blast our zotus
  blastn \
    -db "${db_name}" \
    -outfmt "6 qseqid sseqid staxid ssciname scomname sskingdom pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" \
    ${blast_opt_str} \
    -query ${zotus_fasta} -num_threads ${task.cpus} \
    > blast_result.tsv
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
  blastn -db ${key}_zotus \
    -outfmt "6 qseqid sseqid pident" \
    -out match_list.txt -qcov_hsp_perc 80 \
    -perc_identity 84 -query ${zotus_fasta} \
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
    path 'settings.txt'

  script:
  """
  echo "minimum ratio: ${params.luluMinRatio}" > settings.txt
  echo "minimum ratio type: ${params.luluMinRatioType}" >> settings.txt
  echo "minimum match: ${params.luluMinMatch}" >> settings.txt
  echo "minimum RC: ${params.luluMinRc}" >> settings.txt
  lulu.R \
    -m ${params.luluMinRatio} \
    -t ${params.luluMinRatioType} \
    -a ${params.luluMinMatch} \
    -r ${params.luluMinRc} \
    ${zotu_table} ${match_list} "lulu_zotu_table.tsv" "lulu_zotu_map.tsv" "lulu_result_object.rds"
  """
}

// assign/collapse taxonomy using the R port of the original python script
process collapse_taxonomy {
  label 'r'
  label 'process_single'
  label 'process_more_memory'

  publishDir {
    def td = params.standaloneTaxonomy ? 'standalone_taxonomy' : 'taxonomy'
    "${params.outDir}/${td}/lca/qcov${params.lcaQcov}_pid${params.lcaPid}_diff${params.lcaDiff}"
  }, mode: params.publishMode

  input:
    tuple path(blast_result), path(lineage), path('*')

  output:
    path("lca_taxonomy.tsv"), emit: taxonomy
    path("lca_intermediate.tsv")
    path 'lca_settings.txt'


  script:
  def pf = []
  params.lcaFilterMaxQcov && pf << "--filter-max-qcov"
  params.lcaCaseInsensitive && pf << "--case-insensitive"
  """
  # save settings
  echo "Minimum query coverage %: ${params.lcaQcov}" > lca_settings.txt
  echo "Minimum percent identity: ${params.lcaPid}" >> lca_settings.txt
  echo "Minium percent identity difference: ${params.lcaDiff}" >> lca_settings.txt
  echo "Filter to maximum query coverage: ${params.lcaFilterMaxQcov ? 'yes' : 'no'}" >> lca_settings.txt
  echo "Filter taxa by regex: ${params.lcaTaxonFilter}" >> lca_settings.txt
  echo "Taxon filter case sensitive: ${!params.lcaCaseInsensitive ? 'yes' : 'no'}" >> lca_settings.txt

  collapse_taxonomy.R \
    --qcov ${params.lcaQcov} \
    --pid ${params.lcaPid} \
    --diff ${params.lcaDiff} \
    --merged merged.dmp \
    --nodes nodes.dmp \
    --taxid-lineage taxidlineage.dmp \
    --output lca_taxonomy.tsv \
    --dropped "${params.dropped}" \
    --intermediate "lca_intermediate.tsv" \
    ${pf.join(" ")} \
    --taxon-filter "${params.lcaTaxonFilter}" \
    ${blast_result} ${lineage}
  """
}

// run insect classifier model
process insect {
  label 'insect'
  label 'all_cpus'

  publishDir {
    def offs = String.format("%d",(Integer)num(params.insectOffset))
    def thresh = String.format("%.2f",num(params.insectThreshold))
    def minc = String.format("%d",(Integer)num(params.insectMinCount))
    def ping = String.format("%.2f",num(params.insectPing))
    return "${params.outDir}/taxonomy/insect/thresh${thresh}_offset${offs}_mincount${minc}_ping${ping}"
  }, mode: params.publishMode

  input:
    tuple path(classifier), path('*'), path(zotus)

  output:
    path('insect_taxonomy.tsv'), emit: taxonomy
    path('insect_model.rds')
    path('insect_settings.txt')

  script:
  def offs = String.format("%d",(Integer)num(params.insectOffset))
  def thresh = String.format("%.2f",num(params.insectThreshold))
  def minc = String.format("%d",(Integer)num(params.insectMinCount))
  def ping = String.format("%.2f",num(params.insectPing))

  """
  # record insect settings
  echo "Offset: ${offs}" > insect_settings.txt
  echo "Threshold: ${thresh}" >> insect_settings.txt
  echo "Minimum count: ${minc}" >> insect_settings.txt
  echo "Ping: ${ping}" >> insect_settings.txt

  if [ "${classifier}" != "insect_model.rds" ]; then
    mv ${classifier} insect_model.rds
  fi
  insect.R \
     --cores ${task.cpus} \
     --threshold ${params.insectThreshold} \
     --offset ${params.insectOffset} \
     --min-count ${params.insectMinCount} \
     --ping ${params.insectPing} \
     --lineage rankedlineage.dmp \
     --merged merged.dmp \
     ${zotus} insect_model.rds "insect_taxonomy.tsv"
  """
}

// dummy process to generate published file
process save_taxdump {
  label 'shell'
  label 'process_single'

  publishDir 'output/taxonomy/ncbi_taxdump'

  input:
    path(taxdump)

  output:
    path(taxdump)

  script:
  """
  echo "linking taxdump to publish dir"
  """
}

// use tab-separated sample map to remap sample IDs
process remap_samples {
  label 'shell'
  label 'process_single'

  input:
    tuple val(id), path(reads), path(sample_map)
  output:
    tuple env(new_id), path(reads)

  script:
  if (reads instanceof Path) {
    """
    new_id=\$(awk -F \$'\\t' '\$2 == "${reads}" {print \$1}' ${sample_map})
    if [ -z "\$new_id" ]; then
      new_id="${id}"
    fi
    """
  } else {
    """
    new_id=\$(awk -F \$'\\t' '\$2 == "${reads[0]}" && \$3 == "${reads[1]}" {print \$1}' ${sample_map})
    if [ -z "\$new_id" ]; then
      new_id="${id}"
    fi
    """
  }
}

// un-gzip gzipped files
process unzip {
  label 'shell'
  label 'process_single'

  input:
    tuple val(id), path(reads)
  output:
    tuple val(id), path("*.fastq")

  script:
  """
  z="${reads}"
  # change extension to .fastq regardless of what it currently is
  fn=\$(basename \$z .gz) # get rid of .gz
  fn="\${fn%.*}.fastq"    # get rid of next extension and add .fastq
  gunzip -c \$z > "\$fn"
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
  // if (method == "insect") {
  //   otu = "representative"
  //   c = "--tax-columns \"kingdom,phylum,class,order,family,genus,species,taxon,NNtaxon,NNrank\""
  // }
  """
  phyloseq.R \
    --out phyloseq.rds \
    ${opt.join(" ")} \
    "${zotu_table}" "${taxonomy}" "${metadata}"
  """
}

process finalize {
  label 'r'
  label 'process_single'
  label 'process_more_memory'

  publishDir {
    def td = params.standaloneTaxonomy ? 'standalone_taxonomy/final' : 'final'
    "${params.outDir}/${td}"
  }, mode: params.publishMode

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
  finalize.R \
    --filter "${params.taxonFilter}" \
    --remap "${params.taxonRemap}" \
    --insect "${insect_taxonomy}" \
    --lca "${lca_taxonomy}" \
    --dropped "${params.dropped}" \
    --controls "${params.controls}" \
    --control-action "${params.controlAction}" \
    --control-threshold "${params.controlThreshold}" \
    --decontam-method "${params.decontamMethod}" \
    --concentration "${params.dnaConcentration}" \
    --abundance-threshold "${params.abundanceThreshold}" \
    --rarefaction-method "${params.rarefactionMethod}" \
    --permutations "${params.permutations}" \
    --taxon-priority "${params.taxonPriority}" \
    --curated "${curated_zotu_table}" \
    ${opt.join(" ")} \
    ${zotu_table}
  """
}


// we reuse fastqc/multiqc processes at different steps so they're
// included from an external module
include { fastqc as first_fastqc }    from './modules/modules.nf'
include { fastqc as second_fastqc }   from './modules/modules.nf'
include { multiqc as first_multiqc }  from './modules/modules.nf'
include { multiqc as second_multiqc } from './modules/modules.nf'
include { get_web as get_model } from './modules/modules.nf'
include { get_web as get_taxdb } from './modules/modules.nf'
include { get_web as get_taxdump } from './modules/modules.nf'
include { extract_zip as extract_lineage } from './modules/modules.nf'
include { extract_zip as extract_ncbi } from './modules/modules.nf'
include { extract_targz as extract_taxdb } from './modules/modules.nf'

workflow {
  // make sure our arguments are all in order
  check_params()

  def directions = []

  // do standalone taxonomy assignment
  if (params.standaloneTaxonomy) {
    // build input channels and get appropriate process
    blast_result = Channel.fromPath(params.blastFile, checkIfExists: true)
    // load lineage and (optionally) other taxonomy dump files
    if (!helper.file_exists(params.lcaLineage)) {
      if (!helper.file_exists(params.taxdump)) {
        if (params.taxdump != "") {
          println(colors.yellow("Taxonomy dump archive '${params.taxdump}' does not exist and will be downloaded"))
        }

        Channel.of('https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip') |
          combine(Channel.fromPath('new_taxdump.zip')) |
          get_taxdump |
          set { taxdump }

        taxdump |
          combine(Channel.of('rankedlineage.dmp')) |
          extract_lineage | collect |
          set { lineage }
        taxdump |
          combine(Channel.of('merged.dmp','nodes.dmp','taxidlineage.dmp')) |
          extract_ncbi | collect | toList |
          set { ncbi_dumps }

        taxdump |
          save_taxdump
      } else {
        zip = Channel.fromPath(params.taxdump)
        zip |
          combine(Channel.of('rankedlineage.dmp')) |
          extract_lineage | collect |
          set { lineage }
        zip |
          combine(Channel.of('merged.dmp','nodes.dmp','taxidlineage.dmp')) |
          extract_ncbi | collect | toList |
          set { ncbi_dumps }
      }
    } else {
      // if user has passed a custom ranked lineage dump, use it instead
      lineage = Channel.fromPath(params.lcaLineage)
      ncbi_dumps = Channel.fromPath('nodumps')
    }


    // run taxonomy process
    blast_result |
      combine(lineage) |
      combine(ncbi_dumps) |
      collapse_taxonomy

    // pull out lca table
    collapse_taxonomy.out.taxonomy |
      set { lca_taxonomy }

    // do this part if the zotu table exists
    if (helper.file_exists(params.zotuTable)) {
      zotu_table = Channel.fromPath(params.zotuTable, checkIfExists: true)
      curated_zotu_table = Channel.fromPath("nofile-curated-zotu-table", checkIfExists: false)
      insect_taxonomy = Channel.fromPath("nofile-insect-taxonomy", checkIfExists: false)
      // run it through finalize
      zotu_table |
        combine(curated_zotu_table) |
        combine(lca_taxonomy) |
        combine(insect_taxonomy) |
        finalize
    }
  } else {
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
        // and we try to pull off something from the beginning to use as a sample ID.
        // we also check to see if any of the files end with .gz and mark them as such if they are
        Channel.fromPath(reads_files, checkIfExists: true) |
          map { [it.baseName.tokenize('_')[0],it] } |
          branch {
            gz: it[1] =~ /\.gz$/
            regular: true
          } |
          set { reads }
      } else if (params.paired) {
        // if fwd and rev point to files that exists, just load them directly
        if ( helper.file_exists(params.fwd) && helper.file_exists(params.rev) ) {
          directions = [params.fwd,params.rev]
          Channel.of(params.project) |
            combine(Channel.fromPath(params.fwd,checkIfExists: true)) |
            combine(Channel.fromPath(params.rev,checkIfExists: true)) |
            map { a,b,c -> [a,[b,c]] } |
            branch {
              gz: it[1][0] =~ /\.gz$/ && it[1][1] =~ /\.gz$/
              regular: true
            } |
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
            if (fwd_path[fwd_path.size()-1] != '/') fwd_path += '/'
            if (rev_path[rev_path.size()-1] != '/') rev_path += '/'

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

          // Also, try to enforce the order of r1/r2 reads because fromFilePairs (really, the system)
          // returns the files in alphabetical order, so that if you have a glob like
          // '/reads/{forwards,backwards}*.fastq', the initial result will look like
          // '[ backwards1.fastq, forwards1.fastq ], [ backwards2.fastq, forwards2.fastq ]'
          // Since we've ended up with a glob here, we try to parse the forward/reverse directions
          // from the {f,r} section of the glob
          directions = parse_directions(pattern)
          if (!directions.size()) {
            exit(1,"Unable to determine read direction patterns")
          }

          Channel.fromFilePairs("${pattern}", checkIfExists: true) |
            map { key,f ->
              // enforce read order and make sure we have a key value
              if (key == "") key = params.project
              [key.replaceAll("-","_"),f.sort{a,b -> a =~ /${directions[0]}/ ? -1 : 1}]
            } |
            ifEmpty {
              // bail if we didn't find anything
              exit(1,"No paired reads matched by pattern '${pattern}'. Check command-line options.")
            } |
            branch {
              // separate gzipped reads for unzipping later
              gz: it[1][0] =~ /\.gz$/ && it[1][1] =~ /\.gz$/
              regular: true
            } |
            set { reads }
        }
      } else {
        println(colors.red("Somehow neither ") + colors.bred("--single") + colors.red(" nor ") + colors.bred("--paired") + colors.red(" were passed and we got to this point"))
        println(colors.red("That should not have happened"))
        exit(1)
      }

      // here we decompress any gzipped reads and concatenate
      // them back together with any that weren't gzipped in the first place
      reads.gz |
        // flatten paired-end reads to maximize parallelism in the unzip process
        // [ id, [f1, r2] ] -> [id, r1], [id, r2]
        ( params.paired ? transpose : map { it } ) |
        unzip |
        // regroup paired-end reads back to [ id, [r1, r2] ]
        // and enforce the read direction using groupTuple's sort parameter
        ( params.paired ? groupTuple(sort: {a,b -> a =~ /${directions[0]}/ ? -1 : 1}) : map { it } ) |
        concat ( reads.regular ) |
        set { reads }

      // remap sample IDs if a sample map is provided
      if (params.sampleMap != "") {
        reads |
          combine( Channel.fromPath(params.sampleMap,checkIfExists: true) ) |
          remap_samples |
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
          filter_merge |
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
            ngsfilter |
            set { rfm_barcodes }
        }

        // continue length filtering and whatnot
        rfm_barcodes |
          filter_length |
          map { [it[0], it[1]]} |
          // group together different barcodes because they're
          // concatenated in relabel
          groupTuple |
          // relabel to fasta
          relabel |
          set { relabeled }

        relabeled.result |
          toList | merge_relabeled |
          set { to_dereplicate }

      } else {
        // here, reads are demultiplexed by barcodes, so they're either
        // all in one or two fastq files (depending on single vs paired end)
        // or they're pooled such that barcode pairs are reused across index pairs

        // split the input fastqs to increase parallelism, if requested
        if (params.split) {
          if (params.paired) {
            reads |
              // flatten the reads tuple
              map { key, reads -> [key,reads[0],reads[1]] } |
              // split fastq files
              splitFastq(by: params.splitBy, file: true, pe: true) |
              // rearrange reads tuple so it looks like [key, [R1,R2]]
              // and add the split number to the key
              map { key, read1, read2 -> ["${key}." + file(read1.BaseName).Extension, [read1,read2]] } |
              set { reads }
          } else {
            // in single-end mode we can just split directly
            reads |
              splitFastq(by: params.splitBy, file: true) |
              map { key, readfile -> ["${key}." + file(readfile.BaseName).Extension, readfile] } |
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
          filter_merge |
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
          ngsfilter |
          filter_length |
          split_samples |
          // we have to flatten here because we can get results that look like
          // [[sample1,sample2,sample3],[sample1,sample2,sample3]]
          flatten |
          // this collectFile will merge all individual splits with the same name, so the
          // tuple above turns into [sample1,sample2,sample3]
          collectFile |
          // get rid of the '__split__' business in the filenames
          map { [it.baseName.replaceFirst(/^__split__/,""), it] } |
          // group together different barcodes because they're
          // concatenated in relabel
          groupTuple |
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

    // build the channel, run dereplication, and set to a channel we can use again
    Channel.of(params.project) |
      combine(to_dereplicate) |
      dereplicate |
      set { dereplicated }

    dereplicated.result |
      set { dereplicated }

    // run blast query, unless skipped
    if (params.blast) {
      // def only works on its own line
      // possibly related to NF issue #804: https://github.com/nextflow-io/nextflow/issues/804
      def bdb
      // get $FLOW_BLAST environment variable
      bdb = helper.get_env("FLOW_BLAST")

      // make --blast-db value a list, if it's not already
      def blasts = params.blastDb
      if (!helper.is_list(blasts))
        blasts = [blasts]

      // if $FLOW_BLAST was set and we're not ignoring it, add it to the front of the list
      if (bdb != "" && !params.ignoreBlastEnv)
        blasts = blasts.plus(0,bdb)

      // get unique blast dbs
      blasts = blasts.unique(false)

      // wildcard to capture blast database files
      def wildcard = "{.n*,.[0-9]*.n*}"

      // collect list of files within blast databases
      // and group them by blast db names
      blasts.inject(null,{ b,d ->
        !b ?
          channel.of(file(d).Name) | combine(channel.fromPath("${d}${wildcard}")) :
          b | concat(channel.of(file(d).Name) | combine(channel.fromPath("${d}${wildcard}")))
      }) |
        groupTuple |
        set { blastdb }

      if (!helper.file_exists(params.lcaLineage)) {
        // try to find taxdb files in any of the supplied blast databases
        def db = blasts.collect {
          blastr ->
            file(blastr).Parent.list().findAll {
              it =~ /taxdb\.bt[id]/
            }.collect {
              "${file(blastr).Parent}/${it}"
            }
        }.getAt(0)

        // get taxdb files (either download or from command line)
        if (params.blastTaxdb == "nofile-blast-taxdb") {
          if (db.size() > 0) {
            // if we found something and we don't want something else, use what we found
            Channel.fromPath(db) |
              collect |
              toList |
              set { taxdb }
          } else {
            // download and extract taxdb from ncbi website
            Channel.of('https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz') |
              combine(Channel.fromPath('taxdb.tar.gz')) |
              get_taxdb |
              combine(Channel.of('taxdb.btd','taxdb.bti')) |
              extract_taxdb |
              collect | toList |
              set { taxdb }
          }
        } else {
          Channel.fromPath("${params.blastTaxdb}/taxdb.bt*", checkIfExists: true) |
            collect |
            toList |
            set { taxdb }
        }
      } else {
        Channel.value([[file("taxdb.bti"),file("taxdb.btd")]]) |
          set { taxdb }
      }

      // run the blast query
      dereplicated |
        map { sid, uniques, zotus, zotutable -> zotus } |
        combine(blastdb) |
        combine(taxdb) |
        blast

      // format output directory name for merged blast results
      def pid = String.format("%d",(Integer)num(params.percentIdentity ))
      def evalue = String.format("%.3f",num(params.evalue))
      def qcov = String.format("%d",(Integer)num(params.qcov))
      def blast_dir = "${params.outDir}/blast/pid${pid}_eval${evalue}_qcov${qcov}_max${params.maxQueryResults}"

      // since we're now doing blasts separately for each database, combine the results
      // and store it below each indiviudal database result
      blast.out.result |
        collectFile(name: 'blast_result_merged.tsv', storeDir: blast_dir) |
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

    // get NCBI lineage dump if needed
    if ((params.collapseTaxonomy && params.blast && !helper.file_exists(params.lcaLineage)) || params.insect != false) {
      // get the NCBI ranked taxonomic lineage dump
      if (!helper.file_exists(params.taxdump)) {

        Channel.of('https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip') |
          combine(Channel.fromPath('new_taxdump.zip')) |
          get_taxdump |
          set { taxdump }
        taxdump |
          combine(Channel.of('rankedlineage.dmp')) |
          extract_lineage | collect |
          set { lineage }
        taxdump |
          combine(Channel.of('merged.dmp','nodes.dmp','taxidlineage.dmp')) |
          extract_ncbi | collect | toList |
          set { ncbi_dumps }

        taxdump |
          save_taxdump
      } else {
        zip = Channel.fromPath(params.taxdump)
        zip |
          combine(Channel.of('rankedlineage.dmp')) |
          extract_lineage | collect |
          set { lineage }
        zip |
          combine(Channel.of('merged.dmp','nodes.dmp','taxidlineage.dmp')) |
          extract_ncbi | collect | toList |
          set { ncbi_dumps }
      }
    }

    if (helper.file_exists(params.lcaLineage)) {
      // if user has passed a custom ranked lineage dump, use it instead
      lineage = Channel.fromPath(params.lcaLineage)
      ncbi_dumps = Channel.fromPath('nodumps')

      // if we're doing an insect classification, we still need rankedlineage.dmp from NCBI
      if (params.insect) {
        Channel.of('https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip') |
          combine(Channel.fromPath('new_taxdump.zip')) |
          get_taxdump |
          set { taxdump }
        taxdump |
          combine(Channel.of('rankedlineage.dmp')) |
          extract_lineage | collect |
          set { ncbi_lineage }
      }
    }

    // run the insect classifier, if so desired
    // this should run in parallel with the blast & lulu processes
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
        Channel.of(url) |
          combine(Channel.fromPath('insect_model.rds')) |
          get_model |
          set { classifier }
      }

      // run the insect classification
      classifier |
        // join the correct lineage dump
        ( helper.file_exists(params.lcaLineage) ? combine(ncbi_lineage) : combine(lineage) ) |
        combine(zotus) |
        insect |
        set { insectized }
    }

    // run taxonomy assignment/collapse script if so requested
    if (params.collapseTaxonomy && params.blast) {
      // then we smash it together with the blast results
      // and run the taxonomy assignment/collapser script
      blast_result |
        combine(lineage) |
        combine(ncbi_dumps) |
        collapse_taxonomy
    }

    // prepare for final output
    if (params.collapseTaxonomy && params.blast) {
      collapse_taxonomy.out.taxonomy |
        set { lca_taxonomy }
    } else {
      Channel.fromPath("nofile-lca-taxonomy",checkIfExists: false) |
        set { lca_taxonomy }
    }

    if (params.insect) {
      insect_taxonomy = insectized.taxonomy
    } else {
      insect_taxonomy = Channel.fromPath("nofile-insect-taxonomy", checkIfExists: false)
    }

    if (params.lulu) {
      lulu.out.result |
        map { it[0] } |
        set { curated_zotu_table }
    } else {
      Channel.fromPath("nofile-curated-zotu-table", checkIfExists: false) |
        set { curated_zotu_table }
    }

    if (params.collapseTaxonomy || params.insect) {
      zotu_table |
        combine(curated_zotu_table) |
        combine(lca_taxonomy) |
        combine(insect_taxonomy) |
        finalize
    }

    /* put all the phyloseq stuff down here */
    if (params.phyloseq && helper.file_exists(params.metadata)) {
      def physeq = true
      switch (params.taxonomy) {
        case "lca":
          if (params.collapseTaxonomy) {
            ph_taxonomy = collapse_taxonomy.out.taxonomy
          } else {
            physeq = false
          }
          break
        case "insect":
          if (params.insect) {
            ph_taxonomy = insectized.taxonomy
          } else {
            physeq = false
          }
          break
        case "combined":
          if (params.insect || params.collapseTaxonomy) {
            ph_taxonomy = finalize.out.taxonomy
          } else {
            physeq = false
          }
          break
        default:
          if (helper.file_exists(params.taxonomy)) {
            ph_taxonomy = Channel.fromPath(params.taxonomy, checkIfExists: true)
          } else {
            physeq = false
          }
          break
      }
      if (physeq) {
        metadata = Channel.fromPath(params.metadata, checkIfExists: true)
        seqs = dereplicated.map { sid, uniques, zotus, zotutable -> zotus }
        phyloseq(zotu_table,ph_taxonomy,metadata,seqs)
      }
    }
  }
}
