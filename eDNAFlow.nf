#!/usr/bin/env nextflow24
nextflow.enable.dsl=2

// pull in the helper class
import helper
import colors

/* some global variables */
exec_denoiser = false

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
      --basename ${sample_id}

    mv ${sample_id}.truncated ${sample_id}_trimmed_merged.fastq
    """
  } else {  
    // if reads are paired-end then merge 
    """
    AdapterRemoval --threads ${task.cpus} --file1 ${reads[0]} --file2 ${reads[1]} \
      --collapse --trimns --trimqualities \
      --minquality $params.minQuality \
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
    tuple val(sample_id), val("${barcode.baseName}"), path("*_annotated.fastq") 

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
    tuple val(sample_id), val(barcode_file), path(fastq_file) 
  
  output:
    tuple val(sample_id), path('*_length_filtered.fastq') 

  script:
  // if we're already demultiplexed we probably don't have the forward_tag and reverse_tag annotations
  def p = params.illuminaDemultiplexed ? "" : "-p 'forward_tag is not None and reverse_tag is not None'" 
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
    vsearch --threads 0 --cluster_unoise "${id}_unique.fasta" --centroids "${id}_centroids.fasta" --minsize ${params.minAbundance}	   
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
      usearch -unoise3 "${id}_unique.fasta"  -zotus "${id}_zotus.fasta" -tabbedout "${id}_unique_unoise3.txt" -minsize ${params.minAbundance}
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
      ${params.denoiser} -unoise3 "${id}_unique.fasta"  -zotus "${id}_zotus.fasta" -tabbedout "${id}_unique_unoise3.txt" -minsize ${params.minAbundance}
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

  publishDir "${params.outDir}/blast", mode: params.publishMode 

  input:
    tuple val(sample_id), path(a), path(zotus_fasta), path(zotu_table), path(blast_db), path(custom_db)

  output:
    tuple val(sample_id), path("blast_p${params.percentIdentity}_e${params.evalue}_q${params.qcov}_m${params.maxQueryResults}.tsv"), emit: result

  script:
  def cdb = (String)custom_db
  if (cdb == "NOTHING") {
    cdb = ""
  } else {
    cdb = "${cdb}/${params.customDbName}"
  }
  """

  # this environment variable needs to be there for the taxonomy to show up properly
  export BLASTDB="${blast_db}"
  # blast our zotus
  blastn -task ${params.blastTask} \
    -db "${blast_db}/nt ${cdb}" \
    -outfmt "6 qseqid sseqid staxid ssciname scomname sskingdom pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" \
    -perc_identity ${params.percentIdentity} -evalue ${params.evalue} \
    -best_hit_score_edge 0.05 -best_hit_overhang 0.25 \
    -qcov_hsp_perc ${params.qcov} -max_target_seqs ${params.maxQueryResults} \
    -query ${zotus_fasta} -num_threads ${task.cpus} \
    -out blast_p${params.percentIdentity}_e${params.evalue}_q${params.qcov}_m${params.maxQueryResults}.tsv
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

// retrieve one of the pre-trained insect models from https://github.com/shaunpwilkinson/insect#classifying-sequences
// because of previous sanity checks, we assume the value passed in `model` is a real one
process get_model {
  input:
    val(model)
  output:
    path('insect_model.rds')

  exec:
    m = model.toLowerCase()
    classifier = task.workDir / 'insect_model.rds' 
    url = new URL(helper.insect_classifiers[m])
    url.withInputStream { stream -> classifier << stream }
}

// run insect classifier model
process insect {
  label 'insect'

  publishDir "${params.outDir}/taxonomy/insect", mode: params.publishMode 

  input:
    tuple path(classifier), path(lineage), path(merged), path(zotus), path(zotu_table)

  output:
    tuple path('insect_*.tsv'), path('insect_model.rds'), emit: result

  script:
  """
  if [ "${classifier}" != "insect_model.rds" ]; then
    mv ${classifier} insect_model.rds
  fi
  insect.R \
     --cores ${task.cpus} \
     --threshold ${params.insectThreshold} \
     --offset ${params.insectOffset} \
     --min-count ${params.insectMinCount} \
     --ping ${params.insectPing} \
     --zotu-table ${zotu_table} \
     --lineage ${lineage} \
     --merged ${merged} \
     ${zotus} insect_model.rds "insect_t${params.insectThreshold}_o${params.insectOffset}_m${params.insectMinCount}_p${params.insectPing}.tsv"
  """
}

// retrieve the NCBI ranked taxonomic lineage dump
// also get the info about merged taxids
process get_lineage {
  label 'python3'

  output:
    tuple path('ranked_lineage.tsv'), path('merged.dmp')

  script:
  """
  curl -LO https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
  unzip -p new_taxdump.zip rankedlineage.dmp > ranked_lineage.tsv
  unzip -p new_taxdump.zip merged.dmp > merged.dmp
  rm new_taxdump.zip
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

  // since we changed the way --filter-uncultured works, let's warn the user about it
  if (params.filterUncultured) {
    println(colors.yellow("The parameter --filter-uncultured is no longer recognized. We now filter out uncultured/cloned/whatever sequences by default."))
    println(colors.yellow("To retain those sequences, use the --keep-uncultured argument"))
  }

  // give example of what a demultiplexed FASTA file looks like
  if (params.demuxedExample) {
    helper.demuxed_example()
    exit(0)
  }

  // give a little message that --base-dir is not the way to do it
  if (params.baseDir != "") {
    println(colors.yellow("Parameter ") + colors.byellow("--base-dir") + colors.yellow(" has been replaced by ") + 
            colors.byellow("--reads") + colors.yellow(". Please use ") + colors.byellow("--reads") + colors.yellow(" going forward")) 
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
        if (!params.assignTaxonomy || !params.collapseTaxonomy) {
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
  if (!helper.is_dir(params.blastDb)) {
    println(colors.red("BLAST database must be specified either with the ") + colors.bred("--blast-db") + colors.red(" argument"))
    println(colors.red("or using the \$BLASTDB environment variable. It must point to the directory"))
    println(colors.red("containing the `nt` database (do not include /nt in the path)"))
    exit(1)
  }

  // make sure custom blast database is specified correctly
  if (params.customDbName != "NOTHING" || params.customDb != "NOTHING") {
    if (!helper.is_dir(params.customDb)) {
      println(colors.red("Custom BLAST database must be specified as follows:"))
      println(colors.bred("--custom-db-dir") + colors.red(" <path to custom BLAST db directory>"))
      println(colors.bred("--custom-db") + colors.red(" <name of custom BLAST db (basename of .ndb, etc. files)>"))
      println(colors.red("example:"))
      println(colors.bred("--custom-db-dir") + colors.red(" /storage/blast ") + colors.bred("--custom-db") + colors.red(" test"))
      exit(1)
    }
    if (params.customDbName == "NOTHING") {
      println(colors.red("Name of custom BLAST db must be specified with the ") + colors.bred("--custom-db") + colors.red(" option"))
      exit(1)
    }
  }

  // another blast database sanity check
  if (helper.basename(params.blastDb) == helper.basename(params.customDb)) {
    println(colors.red("Due to the vicissitudes of nextflow internality, the directory names"))
    println(colors.red("of the main and custom BLAST databases must be different."))
    println(colors.red("As specified, both reside in directories called ${helper.basename(params.blastDb)}"))
    exit(1)
  }

  // make sure insect parameter is valid: either a file or one of the pretrained models
  if (params.insect) {
    if (!helper.insect_classifiers.containsKey(params.insect.toLowerCase())) {
      if (!helper.file_exists(params.insect)) {
        println(colors.red("Value passed to ") + colors.bred("--insect") + colors.red(" must be one of the supported builtins or an RDS file"))
        println(colors.red("containing a trained insect classifier model."))
        println(colors.red("See eDNAFlow.nf ") + colors.bred("--help") + colors.red(" for supported builtin models"))
        exit(1)
      }
    }    
  }
}

// we reuse fastqc/multiqc processes at different steps so they're
// included from an external module
include { fastqc as first_fastqc }    from './modules/modules.nf'
include { fastqc as second_fastqc }   from './modules/modules.nf'
include { multiqc as first_multiqc }  from './modules/modules.nf'
include { multiqc as second_multiqc } from './modules/modules.nf'
include { r_lca as collapse_taxonomy_r  } from './modules/modules.nf'
include { r_lca as collapse_taxonomy_lulu_r  } from './modules/modules.nf'
include { py_lca as collapse_taxonomy_py  } from './modules/modules.nf'
include { py_lca as collapse_taxonomy_lulu_py  } from './modules/modules.nf'

workflow {
  // make sure our arguments are all in order
  check_params()


  vsearch = (params.vsearch || params.denoiser == 'vsearch')

  // do standalone taxonomy assignment
  if (params.standaloneTaxonomy) {
    // build input channels and get appropriate process
    tax_process = params.oldTaxonomy ? collapse_taxonomy_py : collapse_taxonomy_r
    zotu_table = Channel.fromPath(params.zotuTable, checkIfExists: true)
    blast_result = Channel.fromPath(params.blastFile, checkIfExists: true)

    // load the ranked lineage and merged channels
    if (!helper.file_exists(params.lineage)) {
      if (params.lineage != "") {
        println(colors.yellow("Lineage file '${params.lineage}' does not exist and will be downloaded"))
      }
      get_lineage |
        set{ lineage }
    } else {
      lineage = Channel.fromPath(params.lineage, checkIfExists: true)
      merged = Channel.fromPath(params.merged, checkIfExists: false)
      lineage | 
        combine(merged) | 
        set { lineage }
    }


    // run taxonomy process
    zotu_table |
      combine(blast_result) |
      combine(lineage) | 
      combine(Channel.of(false)) | 
      tax_process

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
          // fallback to legacy baseDir param if someone forgot
          reads = params.baseDir != "" ? params.baseDir : params.reads
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
      barcodes = Channel.fromPath(params.barcode, checkIfExists: true)

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
        relabel = vsearch ? relabel_vsearch : relabel_usearch

        // run the rest of the pipeline, including the primer mismatch check,
        // length filtering, and smashing together into one file
        reads_filtered_merged | 
          combine(barcodes) | 
          ngsfilter | 
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
        relabel = vsearch ? relabel_vsearch : relabel_usearch

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
      (vsearch ? derep_vsearch : derep_usearch) | set { dereplicated }

      dereplicated.result | 
        set { dereplicated }

    // make blast databases path channels so singularity will automount it properly
    blast_db = Channel.fromPath(params.blastDb, type: 'dir')
    custom_db = Channel.fromPath(params.customDb, type: 'dir', checkIfExists: params.customDb != "NOTHING")

    // run blast query, unless skipped
    if (!params.skipBlast) {
      dereplicated | 
        combine(blast_db) | 
        combine(custom_db) | 
        blast 
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
    if (((params.assignTaxonomy || params.collapseTaxonomy) && !params.skipBlast) || params.insect) {
      // get the NCBI ranked taxonomic lineage dump
      if (!helper.file_exists(params.lineage)) {
        get_lineage |
          set { lineage }
      } else {
        lineage = Channel.fromPath(params.lineage, checkIfExists: true)
        merged = Channel.fromPath(params.merged, checkIfExists: false)
        lineage | 
          combine(merged) | 
          set { lineage }
      }
    }

    // run the insect classifier, if so desired
    // this should run in parallel with the blast & lulu processes
    if (params.insect) {
      // dereplicate returns a tuple, but we only need the zotus fasta
      dereplicated | 
        map { sid, uniques, zotus, zotutable -> [zotus,zotutable] } | 
        set { zotus }
      
      // load the classifier model
      if (helper.file_exists(params.insect)) {
        classifier = Channel.fromPath(params.insect)
      } else {
        // download the classifier model if it's one of the supported ones
        Channel.of(params.insect) | 
          get_model | 
          set { classifier }
      }

      // run the insect classification
      classifier | 
        combine(lineage) | 
        combine(zotus) | 
        insect | set { insectized }
    }

    // run taxonomy assignment/collapse script if so requested
    if ((params.assignTaxonomy || params.collapseTaxonomy) && !params.skipBlast) {
      // here we grab the blast result
      blast.out.result | 
        map { sid, blast_result -> blast_result } | 
        set { blast_result }

    // get the appropriate taxonomy process
      tax_process = params.oldTaxonomy ? collapse_taxonomy_py : collapse_taxonomy_r

      // then we smash it together with the blast results 
      // and run the taxonomy assignment/collapser script
      zotu_table |
        combine(blast_result) |
        combine(lineage) | 
        combine(Channel.of(false)) | 
        tax_process |
        set { taxonomized }
      

      if (!params.skipLulu) {
        // run the taxonomy assignment for lulu-curated zotus
        tax_process_lulu = params.oldTaxonomy ? collapse_taxonomy_lulu_py : collapse_taxonomy_lulu_r
        lulu.out.result | 
          map { zotutable, zotu_map, result_object -> zotutable } | 
          combine(blast_result) |
          combine(lineage) | 
          combine(Channel.of(true)) | 
          tax_process_lulu |
          set { taxonomized_lulu }
      }
    }

    /* put all the phyloseq stuff down here */
    if (params.phyloseq && helper.file_exists(params.metadata)) {
      switch (params.taxonomy) {
        case "lca":
          if (params.assignTaxonomy || params.collapseTaxonomy) {
            taxonomized | 
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
            /* input: tuple path(tax_table), path(zotu_table), path(fasta), path(metadata), val(method) */
            /* tuple val(sample_id), path("${sample_id}_unique.fasta"), path("${sample_id}_zotus.fasta"), path("zotu_table.tsv")  */
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
