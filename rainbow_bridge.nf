#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// pull in the helper class
import helper
import colors

/* some global variables */
exec_denoiser = false


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

  if (params.oldTaxonomy) {
    println(colors.yellow("The parameter --old-taxonomy no longer does anything, since the original python script has been phased out."))
  }

  // give example of what a demultiplexed FASTA file looks like
  if (params.demuxedExample) {
    helper.demuxed_example()
    exit(0)
  }

  if (params.split && params.illuminaDemultiplexed) {
    println(colors.red("Parameters") + colors.bred(" --split ") + colors.red("and") + colors.bred(" --illumina-demultiplexed ") + colors.red("are mutually exclusive"))
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
        if (!params.assignTaxonomy && !params.collapseTaxonomy) {
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
    if (!helper.file_exists(params.zotuTable)) {
      println(colors.red("The supplied zOTU table \"${params.zotuTable}\" does not exist"))  
      exit(1)
    }
  }

  if (params.sampleMap != "" && !helper.file_exists(params.sampleMap)) {
    println(colors.red("The supplied sample map file ${params.sampleMap} does not exist"))
    exit(1)
  }

  // the --vsearch param is no longer supported, but try to catch if someone still uses it
  if (params.containsKey('vsearch')) {
    println(colors.yellow("FYI: vsearch is now used as the default denoiser and the ") + colors.byellow("--vsearch") + colors.yellow(" option is ignored.\nIf you want to use usearch, pass the --usearch option."))
  }

  // if denoiser is an executable, treat it as such
  // otherwise check to make sure it's a valid input
  if (helper.executable(params.denoiser)) {
    exec_denoiser = true
  } else {
    if (!(params.denoiser in ['usearch','usearch32','vsearch'])) {
      println(colors.bred("--denoiser") + colors.red(" must be either 'usearch', 'usearch32', 'vsearch', or a path to an executable (e.g., /opt/sw/bin/usearch64)"))
      exit(1)
    }
  }

  // sanity check, blast database
  if (!params.skipBlast) {

    // if --blast-taxdb is passed, check that it's a directory and that files exist
    if (params.blastTaxdb != "+++__+++") {
      ok = false
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
    bdb = helper.get_env("FLOW_BLAST")

    // make --blast-db param into a list, if it isn't
    blasts = params.blastDb
    if (!helper.is_list(blasts)) 
      blasts = [blasts]

    // if $FLOW_BLAST was set and we're not ignoring it, add it to the front of the list
    if (bdb != "" && !params.ignoreBlastEnv) 
      blasts = blasts.plus(0,bdb)
    
    // get unique vals
    blasts = blasts.unique(false)

    // make sure we've got at least one db
    if (!blasts.size()) {
      println(colors.red("Unless you want to skip the BLAST query with --skip-blast, you must pass at least one value to --blast-db"))
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
  label 'adapterRemoval'
  label 'demux_cpus'

  publishDir "${params.preDir}/trim_merge", mode: params.publishMode

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path('*_trimmed_merged.fastq')

  script:
  if( reads instanceof Path ) {   
    // single end
    """
    AdapterRemoval --threads ${task.cpus} --file1 ${reads} \
      --trimns --trimqualities \
      --minquality ${params.minQuality} \
      --qualitymax ${params.maxQuality} \
      --basename ${sample_id}

    mv ${sample_id}.truncated ${sample_id}_trimmed_merged.fastq
    """
  } else {  
    // if reads are paired-end then merge 
    """
    AdapterRemoval --threads ${task.cpus} --file1 ${reads[0]} --file2 ${reads[1]} \
      --collapse --trimns --trimqualities \
      --minquality $params.minQuality \
      --qualitymax ${params.maxQuality} \
      --minalignmentlength ${params.minAlignLen} \
      --basename ${sample_id}

    mv ${sample_id}.collapsed ${sample_id}_trimmed_merged.fastq  
    """
  }
}

process filter_ambiguous_indices {
  label 'obitools'

  publishDir "${params.preDir}/index_filtered", mode: params.publishMode

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("*_valid_index.fastq") 

  script:
  """
  obigrep --uppercase -D ':[ACGT]+\\+[ACGT]+\$' ${reads} > "${sample_id}_valid_index.fastq"
  """
}

// primer mismatch & sample assignment
// multiple different barcode files are possible
process ngsfilter {
  label 'obitools'

  publishDir "${params.preDir}/ngsfilter", mode: params.publishMode

  input:
    tuple val(sample_id), path(read), path(barcode) 


  output:
    tuple val(sample_id), path("*_annotated.fastq"), val("${barcode.baseName}") 

  script:
  """
  ngsfilter --uppercase -t ${barcode} -e ${params.primerMismatch} -u "${sample_id}_filter_orphans.fastq" ${read} > "${sample_id}_${barcode.baseName}_annotated.fastq"
  """
}

// combine outputs from (possible) multiple barcode files, filter by length
process filter_length {
  label 'obitools'

  publishDir "${params.preDir}/length_filtered", mode: params.publishMode

  input: 
    tuple val(sample_id), path(fastq_file), val(barcode_file) 
  
  output:
    tuple val(sample_id), path('*_length_filtered.fastq') 

  script:
  // if we're already demultiplexed we probably don't have the forward_tag and reverse_tag annotations
  p = params.illuminaDemultiplexed ? "" : "-p 'forward_tag is not None and reverse_tag is not None'" 
  """
  for fastq in ${fastq_file}; do
    obigrep --uppercase -l ${params.minLen} ${p} "\$fastq" >> "${sample_id}_length_filtered.fastq"
  done
  """
}

// for non-demultiplexed runs, split the annotated reads file by samples
process split_samples {
  label 'obitools'

  publishDir "${params.preDir}/split_samples", mode: params.publishMode

  input:
    tuple val(sample_id), path(fastq) 
  
  output:
    path("__split__*.fastq")

  script:
  """
  obisplit --uppercase -p "__split__" -t sample -u "${sample_id}_split_orphans.fastq" $fastq
  """
}

// relabel files for vsearch
process relabel_vsearch {
  label 'vsearch'

  publishDir "${params.preDir}/relabeled", mode: params.publishMode

  input:
    tuple val(sample_id), path(fastq)
  output:
    path('*_relabeled.fasta'), emit: result
    path 'settings.txt'


  script:
  """
  echo "denoiser: vsearch" > settings.txt
  # this may or may not be necessary anymore, but it seems like a good sanity check
  # since this will fail on empty files
  if [ -s "${fastq}" ]; then 
    vsearch --threads 0 --fastx_filter ${fastq} --relabel "${sample_id}." --fastaout - | \
      awk '/^>/ {print;} !/^>/ {print(toupper(\$0))}' > "${sample_id}_relabeled.fasta"
  else
    touch "${sample_id}_relabeled.fasta"
  fi
  """
}

// relabel files for usearch
process relabel_usearch {
  label 'usearch'

  publishDir "${params.preDir}/relabeled", mode: params.publishMode

  input:
    tuple val(sample_id), path(fastq)
  output:
    path('*_relabeled.fasta'), emit: result
    path 'settings.txt'


  script:
  if (!exec_denoiser) {
    """
    echo "denoiser: usearch" > settings.txt
    if [ -s "${fastq}" ]; then 
      # we have to convert everything to uppercase because obisplit --uppercase is broken
      usearch -fastx_relabel ${fastq} -prefix "${sample_id}." -fastaout /dev/stdout | tail -n+7 | \
        awk '/^>/ {print;} !/^>/ {print(toupper(\$0))}' > "${sample_id}"_relabeled.fasta 
    else
      touch "${sample_id}_relabeled.fasta"
    fi
    """
  } else if (exec_denoiser) {
    """
    echo "denoiser: ${params.denoiser}" > settings.txt
    for files in ${fastqs}; do
      label=\$(echo \$files | cut -d '/' -f 3 | cut -d '.' -f 1)
      ${params.denoiser} -fastx_relabel \$files -prefix \${label}. -fastqout \${label}.relabeled.fastq
    done

    for files in *.relabeled.fastq; do
      name=\$(echo \$files | cut -d '/' -f '2' | cut -d '.' -f 1)
      echo \${name} >> CountOfSeq.txt
      grep "^@\${name}" \$files | wc -l >> CountOfSeq.txt
    done

    cat *.relabeled.fastq > "${sample_id}_QF_Dmux_minLF_relabeled4Usearch.fastq"

    ${params.denoiser} -fastx_get_sample_names *_relabeled4Usearch.fastq -output sample.txt

    ${params.denoiser} -fastq_filter *_relabeled4Usearch.fastq -fastaout ${sample_id}_upper.fasta

    # awk '/^>/ {print(\$0)}; /^[^>]/ {print(toupper(\$0))}' *.fasta > ${sample_id}_upper.fasta
    """
  } else {
    """
    echo "we were passed a mode that wasn't usearch, usearch64, or vsearch"
    exit 1
    """
  }
}

// dereplication, zOTUs creation, zOTU table creation (vsearch version)
process derep_vsearch {
  label 'vsearch'
  label 'all_cpus'

  publishDir "${params.outDir}/zotus", mode: params.publishMode

  input:
    tuple val(id), path(relabeled_merged)

  output:
    tuple val(id), path("${id}_unique.fasta"), path("${id}_zotus.fasta"), path("zotu_table.tsv"), emit: result
    path 'settings.txt'

  script:
  """
  echo "denoiser: vsearch" > settings.txt
  echo "minimum sequence abundance: ${params.minAbundance}" >> settings.txt
  # pass 0 threads to use all available cores
  # steps:
  # 1. get unique sequence variants
  # 2. run denoising algorithm
  # 3. get rid of chimeras
  # 4. match original sequences to zotus by 97% identity
  if [ -s "${relabeled_merged}" ]; then 
    vsearch --threads 0 --derep_fulllength ${relabeled_merged} --sizeout --output "${id}_unique.fasta"
    vsearch --threads 0 --cluster_unoise "${id}_unique.fasta" --centroids "${id}_centroids.fasta" --minsize ${params.minAbundance} --unoise_alpha ${params.alpha}
    vsearch --threads 0 --uchime3_denovo "${id}_centroids.fasta" --nonchimeras "${id}_zotus.fasta" --relabel Zotu 
    vsearch --threads 0 --usearch_global ${relabeled_merged} --db "${id}_zotus.fasta" --id 0.97 --otutabout zotu_table.tsv
  else
    >&2 echo "Merged FASTA is empty. Did your PCR primers match anything?"  
    exit 1
  fi
  """
}

// dereplication, etc. using usearch
process derep_usearch {
  label 'usearch'
  label 'all_cpus'

  publishDir "${params.outDir}/zotus", mode: params.publishMode

  input:
    tuple val(id), path(relabeled_merged) 

  output:
    tuple val(id), path("${id}_unique.fasta"), path("${id}_zotus.fasta"), path("zotu_table.tsv"), emit: result
    path 'settings.txt'

  script:
  if (!exec_denoiser)
  {
    """
    echo "denoiser: usearch" > settings.txt
    echo "minimum sequence abundance: ${params.minAbundance}" >> settings.txt
    # steps:
    # 1. get unique sequences
    # 2. run denoising & chimera removal
    # 3. generate zotu table
    if [ -s "${relabeled_merged}" ]; then
      usearch -fastx_uniques ${relabeled_merged} -sizeout -fastaout "${id}_unique.fasta"
      usearch -unoise3 "${id}_unique.fasta"  -zotus "${id}_zotus.fasta" -tabbedout "${id}_unique_unoise3.txt" -minsize ${params.minAbundance} -unoise_alpha ${params.alpha}
      usearch -otutab ${relabeled_merged} -zotus ${id}_zotus.fasta -otutabout zotu_table.tsv -mapout zmap.txt
    else
      >&2 echo "${colors.bred('Merged FASTA is empty. Did your PCR primers match anything?')}"  
      exit 1
    fi
    """
  } else if (exec_denoiser) {
    """
    echo "denoiser: ${params.denoiser}" > settings.txt
    echo "minimum sequence abundance: ${params.minAbundance}" >> settings.txt
    if [ -s "${relabeled_merged}" ]; then
      ${params.denoiser} -fastx_uniques ${relabeled_merged} -sizeout -fastaout "${id}_unique.fasta"
      ${params.denoiser} -unoise3 "${id}_unique.fasta"  -zotus "${id}_zotus.fasta" -tabbedout "${id}_unique_unoise3.txt" -minsize ${params.minAbundance} -unoise_alpha ${params.alpha}
      ${params.denoiser} -otutab ${relabeled_merged} -zotus ${id}_zotus.fasta -otutabout zotu_table.tsv -mapout zmap.txt
    else
      >&2 echo "${colors.bred('Merged FASTA is empty. Did your PCR primers match anything?')}"  
      exit 1
    fi
    """
  } else {
    """
    echo "we were passed a mode that wasn't usearch, usearch64, or vsearch"
    exit 1
    """
  }
}

// run blast query
process blast {
  label 'blast'

  publishDir { 
    pid = String.format("%d",(Integer)num(params.percentIdentity ))
    evalue = String.format("%.3f",num(params.evalue))
    qcov = String.format("%d",(Integer)num(params.qcov))
    return "${params.outDir}/blast/pid${pid}_eval${evalue}_qcov${qcov}_max${params.maxQueryResults}/${db_name}" 
  }, mode: params.publishMode 

  input:
    tuple path(zotus_fasta), val(db_name), path(db_files), path(taxdb)

  output:
    path("blast_result.tsv"), emit: result
    path 'blast_settings.txt'

  script:

  // format settings values
  pid = String.format("%d",(Integer)num(params.percentIdentity))
  evalue = String.format("%.3f",num(params.evalue))
  qcov = String.format("%d",(Integer)num(params.qcov))           
  """
  # record blast settings
  echo "Percent identity: ${pid}" > blast_settings.txt
  echo "e-value: ${evalue}" >> blast_settings.txt
  echo "Query qoverage: ${qcov}" >> blast_settings.txt
  echo "Max. target sequences: ${params.maxQueryResults}" >> blast_settings.txt

  # set BLASTDB to local working directory
  export BLASTDB=.

  # blast our zotus
  blastn -task ${params.blastTask} \
    -db "${db_name}" \
    -outfmt "6 qseqid sseqid staxid ssciname scomname sskingdom pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" \
    -perc_identity ${params.percentIdentity} -evalue ${params.evalue} \
    -best_hit_score_edge 0.05 -best_hit_overhang 0.25 \
    -qcov_hsp_perc ${params.qcov} -max_target_seqs ${params.maxQueryResults} \
    -query ${zotus_fasta} -num_threads ${task.cpus} \
    -out blast_result.tsv
  """
}

// make custom blast database for LULU curation
process lulu_blast {
  label 'blast'

  input:
    tuple val(sample_id), path(zotus_fasta), path(zotu_table)
  
  output:
    tuple val(sample_id), path('match_list.txt'), path(zotu_table)

  script:
  """
  # blast zotus against themselves to create the match list LULU needs
  makeblastdb -in ${zotus_fasta} -parse_seqids -dbtype nucl -out ${sample_id}_zotus
  blastn -db ${sample_id}_zotus \
    -outfmt "6 qseqid sseqid pident" \
    -out match_list.txt -qcov_hsp_perc 80 \
    -perc_identity 84 -query ${zotus_fasta} \
    -num_threads ${task.cpus}
  """
}

// LULU curation
process lulu {
  label 'r'

  publishDir "${params.outDir}/lulu", mode: params.publishMode

  input:
    tuple val(sample_id), path(match_list), path(zotu_table)

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


// run insect classifier model
process insect {
  label 'insect'

  publishDir { 
    offs = String.format("%d",(Integer)num(params.insectOffset))
    thresh = String.format("%.2f",num(params.insectThreshold))
    minc = String.format("%d",(Integer)num(params.insectMinCount))
    ping = String.format("%.2f",num(params.insectPing))
    return "${params.outDir}/taxonomy/insect/thresh${thresh}_offset${offs}_mincount${minc}_ping${ping}"
  }, mode: params.publishMode 

  input:
    tuple path(classifier), path(lineage), path(merged), path(zotus)

  output:
    tuple path('insect_taxonomy.tsv'), path('insect_model.rds'), emit: result
    path 'insect_settings.txt'

  script:
  offs = String.format("%d",(Integer)num(params.insectOffset))
  thresh = String.format("%.2f",num(params.insectThreshold))
  minc = String.format("%d",(Integer)num(params.insectMinCount))
  ping = String.format("%.2f",num(params.insectPing))

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
     --lineage ${lineage} \
     --merged ${merged} \
     ${zotus} insect_model.rds "insect_taxonomy.tsv"
  """
}

// extract the NCBI blast taxonomy database
process extract_taxdb {
  label 'python3'

  input:
    tuple path(taxdb), val(to_extract)
  
  output:
    path(to_extract)

  script:
  """
  tar -zxvf taxdb.tar.gz ${to_extract}
  """
}

// extract arbitrary files from a zip archive
process extract_taxonomy {
  label 'python3'

  input:
    tuple path(zipfile), val(f)
  output:
    path(f)
  
  script:
  """
  unzip -p ${zipfile} ${f} > ${f}
  """
}

// dummy process to generate published file
process save_taxdump {
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
  label 'unzip'

  input:
    tuple val(id), path(reads)
  output: 
    tuple val(id), path("*.fastq")

  script:
  """
  zips=( ${reads} )
  for z in "\${zips[@]}"; do
    gunzip -c \$z > \${z%.gz}
  done
  """
}

process phyloseq {
  label 'r'

  publishDir "${params.outDir}/phyloseq", mode: params.publishMode

  input:
    tuple path(tax_table), path(zotu_table), path(fasta), path(metadata), val(method)
  
  output:
    path("phyloseq_${method}.rds")

  script:
  otu = "OTU"
  t = params.noTree ? "--no-tree" : ""
  o = params.optimizeTree ? "--optimize" : ""
  c = params.taxColumns != "" ? "--tax-columns ${params.taxColumns}" : ""
  if (method == "insect") {
    otu = "representative"
    c = "--tax-columns \"kingdom,phylum,class,order,family,genus,species,taxon,NNtaxon,NNrank\""
  }
  """
  phyloseq.R \
    --taxonomy "${tax_table}" \
    --otu "${otu}" \
    --otu-table "${zotu_table}" \
    --cores ${task.cpus} \
    --fasta "${fasta}" \
    --metadata "${metadata}" \
    --phyloseq "phyloseq_${method}.rds" ${t} ${o} ${c}
  """
}

process finalize {
  label 'r'
  
  publishDir "${params.outDir}/final", mode: params.publishMode

  input:
    tuple path(zotu_table), path(curated_zotu_table), path(lca_taxonomy), path(insect_taxonomy)

  output:
    tuple path("zotu_table_raw.tsv"), path("taxonomy_raw.tsv"), path("zotu_table_final*.tsv")

  script:
  opt = []
  if (params.abundanceFilter) {
    opt.add("--abundance-filter")
  }
  if (params.rarefy) {
    opt.add("--rarefy")
  }
  if (params.filterMinimum) {
    opt.add("--filter-min")
  }
  """
  finalize.R \
    --filter "${params.taxonFilter}" \
    --remap "${params.taxonRemap}" \
    --insect "${insect_taxonomy}" \
    --lca "${lca_taxonomy}" \
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
include { lca as collapse_taxonomy  } from './modules/modules.nf'
include { get_web as get_model } from './modules/modules.nf'
include { get_web as get_taxdb } from './modules/modules.nf'
include { get_web as get_taxdump } from './modules/modules.nf'

workflow {
  // make sure our arguments are all in order
  check_params()

  lca = params.collapseTaxonomy || params.assignTaxonomy
  usearch = (params.usearch || params.denoiser == 'usearch')

  // do standalone taxonomy assignment
  if (params.standaloneTaxonomy) {
    // build input channels and get appropriate process
    zotu_table = Channel.fromPath(params.zotuTable, checkIfExists: true)
    blast_result = Channel.fromPath(params.blastFile, checkIfExists: true)

    // load the ranked lineage and merged channels
    if (!helper.file_exists(params.taxdump)) {
      if (params.taxdump != "") {
        println(colors.yellow("Taxonomy dump archive '${params.taxdump}' does not exist and will be downloaded"))
      }

      Channel.of('https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip') | 
        combine(Channel.fromPath('new_taxdump.zip')) | 
        get_taxdump | 
        set { taxdump }
      taxdump | 
        combine(Channel.of('rankedlineage.dmp','merged.dmp','nodes.dmp')) |
        extract_taxonomy | collect | toList |
        set { lineage }

      taxdump | 
        save_taxdump
    } else {
      Channel.fromPath(params.taxdump) | 
        combine(Channel.of('rankedlineage.dmp','merged.dmp','nodes.dmp')) |
        extract_taxonomy | collect | toList |
        set { lineage }
    }


    // run taxonomy process
    blast_result |
      combine(lineage) | 
      collapse_taxonomy

  } else {
    if (!helper.file_exists(params.demuxedFasta)) {
      if (params.single) {
        // here we load whatever was passed as the --reads option
        // if it's a glob, we get a list of files. if it's just one, we get just one
        // and we try to pull off something from the beginning to use as a sample ID. 
        // we also check to see if any of the files end with .gz and mark them as such if they are
        // (that's what branch does)
        Channel.fromPath(params.reads, checkIfExists: true) |
          map { [it.baseName.tokenize('_')[0],it] } |
          branch { 
            gz: it[1] =~ /\.gz$/
            regular: true
          } |
          set { reads }
      } else if (params.paired) { 
        // if fwd and rev point to files that exists, just load them directly
        // I guess these can probably be globs if you want them to be, but that
        // might not work properly as currently written
        if ( helper.file_exists(params.fwd) && helper.file_exists(params.rev) ) {
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
          // TODO: if params.reads is *not* a directory, should we just treat it as a glob and pass it directly into
          // fromFilePairs?
          reads = params.reads
          // otherwise make a glob pattern to find where the fwd/rev reads live
          // this may get a little complex, so there are two possible options:
          // files must have some way of stating which direction they are (e.e., R1/R2)
          // fwd/rev reads may optionally be in separate subdirectories but those must both be subdirectories
          // at the same level. So you can have /dir/R1 and /dir/R2, but you CAN'T have /dir1/R1 and /dir2/R2
          if (params.fwd != "" && params.rev != "") {
            pattern = "${reads}/{${params.fwd},${params.rev}}/*{${params.r1},${params.r2}}*.fastq*"
          } else {
            pattern = "${reads}/*{${params.r1},${params.r2}}*.fastq*"
          }

          // load the reads and make sure they're in the right order because the
          // files will be sorted alphabetically. this is an extremely niche issue
          // because file are basically always R1/R2 but what if they're
          // forwards/backwards or something like that?

          // we also replace dashes with underscores in the sample ID because at least
          // vsearch will cut on that character and not treat it as a complete sample id

          // I'm not even certain the order matters, but I think it does because we 
          // send them to --file1 and --file2 of AdapterRemoval
          Channel.fromFilePairs(pattern, checkIfExists: true) | 
            map { key,f -> [key.replaceAll("-","_"),f[0] =~ /${params.r1}/ ? [f[0],f[1]] : [f[1],f[0]]]  } | 
            map { key,f -> key == "" ? [params.project,f] : [key,f] } | 
            branch { 
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

      // here we unzip the zipped files (if any) and merge them back together
      // with the unzipped files (if any)
      reads.gz | 
        unzip |
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
      barcodes = Channel.fromPath(params.barcode, checkIfExists: !params.skipPrimerMatch)

      // if the sequences are already demultiplexed by illumina, we'll
      // process them separately, including optionally attempting to remove ambiguous indices
      // and ultimately smash them together for vsearch/usearch to do the dereplication
      if (params.illuminaDemultiplexed) {

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

        // get relabeling process
        relabel = usearch ? relabel_usearch : relabel_vsearch

        // with or without the primer mismatch check, do the
        // length filtering and smash results together into one file
        reads_filtered_merged |
          combine(barcodes) |
          set { rfm_barcodes }

        // only run ngsfilter if we have primers
        if(!params.skipPrimerMatch) {
          rfm_barcodes |
            ngsfilter |
            set { rfm_barcodes }
        }

        // continue length filtering and whatnot
        rfm_barcodes |
          groupTuple | 
          filter_length |
          relabel |
          set { relabeled }

        relabeled.result | 
          collectFile(name: "${params.project}_relabeled.fasta", storeDir: "${params.preDir}/merged") |
          set { to_dereplicate } 

      } else {
        // this is where reads are all in one fwd or fwd/rev file and have NOT been demultiplexed

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
          if (params.split) {
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

        // do after-filtering fastqc step
        if (params.fastqc) {
          Channel.of("filtered") | 
            combine(reads_filtered_merged) |
            second_fastqc
          // again run multiqc if split
          if (params.split) {
            second_fastqc.out | 
              collect(flat: true) | 
              toList |
              combine(Channel.of("filtered")) | 
              second_multiqc
          }
        }

        // get relabeling process
        relabel = usearch ? relabel_usearch : relabel_vsearch

        // run the rest of the pipeline, including demultiplexing, length filtering,
        // splitting, and recombination for dereplication
        reads_filtered_merged | 
          combine(barcodes) |
          ngsfilter | 
          groupTuple | 
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
          // relabel to fasta
          relabel | 
          set { relabeled }
        // collect to single relabeled fasta
        relabeled.result | 
          collectFile(name: "${params.project}_relabeled.fasta", storeDir: "${params.preDir}/merged") |
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
      (usearch ? derep_usearch : derep_vsearch) | set { dereplicated }

    dereplicated.result | 
      set { dereplicated }

    // run blast query, unless skipped
    if (!params.skipBlast) {
      // get $FLOW_BLAST environment variable
      bdb = helper.get_env("FLOW_BLAST")

      // make --blast-db value a list, if it's not already
      blasts = params.blastDb
      if (!helper.is_list(blasts)) 
        blasts = [blasts]

      // if $FLOW_BLAST was set and we're not ignoring it, add it to the front of the list
      if (bdb != "" && !params.ignoreBlastEnv) 
        blasts = blasts.plus(0,bdb)

      // get unique vals
      blasts = blasts.unique(false)

      // make wildcards for blast database files
      blast_db_files = blasts.collect {
        "${it}*.n*"
      }

      // construct input channels

      // collect list of files within blast databases
      Channel.fromPath(blast_db_files, checkIfExists: true) | 
        map { [ it.Name.toString().replaceAll(/(\.[0-9]+)?\.n..$/,''), it  ]  } | 
        groupTuple | 
        set { blastdb }

      // try to find taxdb files in any of the supplied blast databases
      db = blasts.collect {
        blastr -> 
          file(blastr).Parent.list().findAll {
            it =~ /taxdb\.bt[id]/
          }.collect {
            "${file(blastr).Parent}/${it}"
          }
      }.getAt(0)


      // get taxdb files (either download or from command line)
      if (params.blastTaxdb == "+++__+++") {
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

      // run the blast query
      dereplicated | 
        map { sid, uniques, zotus, zotutable -> zotus } | 
        combine(blastdb) |
        combine(taxdb) | 
        blast 

      // format output directory name for merged blast results
      pid = String.format("%d",(Integer)num(params.percentIdentity ))
      evalue = String.format("%.3f",num(params.evalue))
      qcov = String.format("%d",(Integer)num(params.qcov))
      blast_dir = "${params.outDir}/blast/pid${pid}_eval${evalue}_qcov${qcov}_max${params.maxQueryResults}" 

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
    if (!params.skipLulu) {
      dereplicated | 
        // get zotus and sample id
        map { sid, uniques, zotus, zotutable -> [sid,zotus,zotutable] } | 
        lulu_blast | 
        lulu
    }

    // get NCBI lineage dump if needed
    if ((lca && !params.skipBlast) || params.insect) {
      // get the NCBI ranked taxonomic lineage dump
      if (!helper.file_exists(params.taxdump)) {

        Channel.of('https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip') | 
          combine(Channel.fromPath('new_taxdump.zip')) | 
          get_taxdump | 
          set { taxdump }
        taxdump | 
          combine(Channel.of('rankedlineage.dmp','merged.dmp','nodes.dmp')) |
        extract_taxonomy | collect | toList |
        set { lineage }

        taxdump | 
          save_taxdump
      } else {
        Channel.fromPath(params.taxdump) | 
          combine(Channel.of('rankedlineage.dmp','merged.dmp','nodes.dmp')) |
        extract_taxonomy | collect | toList |
        set { lineage }
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
        m = params.insect.toLowerCase()
        url = helper.insect_classifiers[m]
        Channel.of(url) | 
          combine(Channel.fromPath('insect_model.rds')) |
          get_model | 
          set { classifier }
      }

      // run the insect classification
      classifier | 
        combine(lineage) | 
        combine(zotus) | 
        insect | 
        set { insectized }
    }

    // run taxonomy assignment/collapse script if so requested
    if (lca && !params.skipBlast) {
      // then we smash it together with the blast results 
      // and run the taxonomy assignment/collapser script
      blast_result |
        combine(lineage) | 
        collapse_taxonomy |
        set { taxonomized }
    }

    // prepare for final output
    if (lca && !params.skipBlast) {
      taxonomized.result | 
        map { it[1] } | 
        set { lca_taxonomy }
    } else {
      Channel.fromPath("NOTADANGFILE.nothing.lca",checkIfExists: false) | 
        set { lca_taxonomy } 
    }

    if (params.insect) {
      insectized.result | 
        map { it[0] } |
        set { insect_taxonomy }
    } else {
      Channel.fromPath("NOTADANGFILE.nothing.insect", checkIfExists: false) | 
        set { insect_taxonomy }
    }

    if (!params.skipLulu) {
      lulu.out.result |
        map { it[0] } |
        set { curated_zotu_table }
    } else {
      Channel.fromPath("NOTADANGFILE.nothing.lulu", checkIfExists: false) | 
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
      switch (params.taxonomy) {
        case "lca":
          if (lca) {
            taxonomized.result | 
              map { it[1] } |
              combine(zotu_table) | 
              combine( dereplicated | map { sid, uniques, zotus, zotutable -> zotus } ) |
              combine( Channel.fromPath(params.metadata, checkIfExists: true) ) | 
              combine( Channel.of("lca") ) |
              phyloseq 
          }
          break
        case "insect":
          if (params.insect != false) {
            insectized.result | 
              map { it[0] } |
              combine(zotu_table) | 
              combine( zotus ) |
              combine( Channel.fromPath(params.metadata, checkIfExists: true) ) | 
              combine( Channel.of("insect") ) |
              phyloseq
          }
          break
        default:
          if (helper.file_exists(params.taxonomy)) {
            Channel.fromPath(params.taxonomy, checkIfExists: true) | 
              combine(zotu_table) | 
              combine( dereplicated | map { sid, uniques, zotus, zotutable -> zotus } ) |
              combine( Channel.fromPath(params.metadata, checkIfExists: true) ) | 
              combine( Channel.of("file") ) |
              phyloseq
          } 
          break
      }
    }
  }
}

