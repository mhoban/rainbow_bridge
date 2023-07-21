#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// pull in the helper class
import helper

/* some global variables */
exec_denoiser = false

// trim and (where relevant) merge paired-end reads
// this gets called if the files were already demultiplexed
process filter_merge {
  label 'adapterRemoval'
  label 'demux_cpus'

  publishDir { params.illuminaDemultiplexed ? '01_demultiplexed_filter_merge' : '01_filter_merge' }, mode: params.publishMode

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

    split_id=""
    if ${params.split}; then
      split_id="_${reads.baseName}"
    fi

    mv ${sample_id}.truncated ${sample_id}\${split_id}_trimmed_merged.fastq
    """
  } else {  
    // if reads are paired-end then merge 
    """
    AdapterRemoval --threads ${task.cpus} --file1 ${reads[0]} --file2 ${reads[1]} \
      --collapse --trimns --trimqualities \
      --minquality $params.minQuality \
      --minalignmentlength ${params.minAlignLen} \
      --basename ${sample_id}

    split_id=""
    if ${params.split}; then
      split_id="_${reads[0].baseName}"
    fi

    mv ${sample_id}.collapsed ${sample_id}\${split_id}_trimmed_merged.fastq  
    """
  }
}

process filter_ambiguous_tags {
  label 'obitools'

  publishDir '01a_ambiguous_tags_filtered', mode: params.publishMode

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("*_good_tags.fastq") 

  script:
  """
  obigrep --uppercase -D ':[ACGT]+\\+[ACGT]+\$' ${reads} > "${sample_id}_good_tags.fastq"
  """
}

// annotate sequences that have already been demultiplexed by illumina
// we probably don't actually need this but it's being left in for now
process annotate_demultiplexed {
  label 'obitools'

  publishDir '0x_annotated', mode: params.publishMode

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("*_annotated.fastq")

  script:
  if (params.removeAmbiguousTags) {
    """
    # optionally remove ambiguous tags and annotate sample names
    # this assumes that illumina-demultiplexed sequence files have their tag sequences in the header
    # line looking like [AGCT]+:[AGCT]+
    obiannotate --uppercase -S illumina_sample:"'${sample_id}'" ${reads} | obigrep --uppercase -D ':[ACGT]+\\+[ACGT]+\$' > "${sample_id}_annotated.fastq"
    """
  } else {
    """
    # annotate sample names
    obiannotate --uppercase -S illumina_sample:"'${sample_id}'" ${reads} > "${sample_id}_annotated.fastq"
    """
  }
}

// primer mismatch & sample assignment
// multiple different barcode files are possible
process ngsfilter {
  label 'obitools'

  publishDir "02_ngsfilter"

  input:
    tuple val(sample_id), path(read), path(barcode) 


  output:
    tuple val(sample_id), val("${barcode.baseName}"), path("*_Dmux.fastq") 

  script:
  """
  ngsfilter --uppercase -t ${barcode} -e ${params.primerMismatch} -u "orphan.fastq" ${read} > "${sample_id}_${read.baseName}_${barcode.baseName}_QF_Dmux.fastq"
  """
}

// combine outputs from (possible) multiple barcode files, filter by length
process filter_length {
  label 'obitools'

  publishDir "03_length_filtered"

  input: 
    tuple val(sample_id), val(barcode_file), path(fastq_file) 
  
  output:
    tuple val(sample_id), path('*_QF_Dmux_minLF.fastq') 

  script:
  // if we're already demultiplexed we probably don't have the forward_tag and reverse_tag annotations
  def p = params.illuminaDemultiplexed ? "" : "-p 'forward_tag is not None and reverse_tag is not None'" 
  """
  # cat ${fastq_file} > "${sample_id}_QF_Dmux.fastq" 
  obigrep --uppercase -l ${params.minLen} ${p} "${fastq_file}" > "${sample_id}_${fastq_file.baseName}_QF_Dmux_minLF.fastq"
  """
}

// for non-demultiplexed runs, split the annotated reads file by samples
process split_samples {
  label 'obitools'

  publishDir "04_splitSamples_${sample_id}", mode: params.publishMode

  input:
    tuple val(sample_id), path(fastq) 
  
  output:
    path("__split__*.fastq")

  script:
  // TODO: we probably don't need this anymore since we're not doing this annotation now
  def sample_tag = params.illuminaDemultiplexed ? "illumina_sample" : "sample"
  """
  obisplit --uppercase -p "__split__" -t ${sample_tag} -u "orphans.fastq" $fastq
  """
}

// relabel files for vsearch
process relabel_vsearch {
  label 'vsearch'

  publishDir { params.illuminaDemultiplexed ? "04_relabel_vsearch" : "05_relabel_vsearch" }, mode: params.publishMode

  input:
    tuple val(sample_id), path(fastq)

  output:
    path('*_relabeled.fasta')

  script:
  """
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

  publishDir { params.illuminaDemultiplexed ? "04_relabel_usearch" : "05_relabel_usearch" }, mode: params.publishMode

  input:
    tuple val(sample_id), path(fastq) 

  output:
    path('*_relabeled.fasta')

  script:
  if (!exec_denoiser) {
    """
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

  publishDir { params.illuminaDemultiplexed ? "06_derep_vsearch" : "07_derep_vsearch" }, mode: params.publishMode

  input:
    tuple val(sample_id), path(upper_fasta) 

  output:
    tuple val(sample_id), path("${sample_id}_Unq.fasta"), path("${sample_id}_zotus.fasta"), path("zotuTable.txt") 

  script:
  """
  # pass 0 threads to use all available cores
  # steps:
  # 1. get unique sequence variants
  # 2. run denoising algorithm
  # 3. get rid of chimeras
  # 4. match original sequences to zotus by 97% identity
  vsearch --threads 0 --derep_fulllength ${upper_fasta} --sizeout --output "${sample_id}_Unq.fasta"
  vsearch --threads 0 --cluster_unoise "${sample_id}_Unq.fasta" --centroids "${sample_id}_centroids.fasta" --minsize ${params.minAbundance}	   
  vsearch --threads 0 --uchime3_denovo "${sample_id}_centroids.fasta" --nonchimeras "${sample_id}_zotus.fasta" --relabel Zotu 
  vsearch --threads 0 --usearch_global ${upper_fasta} --db "${sample_id}_zotus.fasta" --id 0.97 --otutabout zotuTable.txt
  """
}

// dereplication, etc. using usearch
process derep_usearch {
  label 'usearch'
  label 'all_cpus'

  publishDir { params.illuminaDemultiplexed ? "06_derep_usearch" : "07_derep_usearch" }, mode: params.publishMode

  input:
    tuple val(sample_id), path(upper_fasta) 

  output:
    tuple val(sample_id), path("${sample_id}_Unq.fasta"), path("${sample_id}_zotus.fasta"), path("zotuTable.txt") 

  script:
  if (!exec_denoiser)
  {
    """
    # steps:
    # 1. get unique sequences
    # 2. run denoising & chimera removal
    # 3. generate zotu table
    usearch -fastx_uniques ${upper_fasta} -sizeout -fastaout "${sample_id}_Unq.fasta"
    usearch -unoise3 "${sample_id}_Unq.fasta"  -zotus "${sample_id}_zotus.fasta" -tabbedout "${sample_id}_Unq_unoise3.txt" -minsize ${params.minAbundance}
    usearch -otutab ${upper_fasta} -zotus ${sample_id}_zotus.fasta -otutabout zotuTable.txt -mapout zmap.txt
    """
  } else if (exec_denoiser) {
    """
    ${params.denoiser} -fastx_uniques ${upper_fasta} -sizeout -fastaout "${sample_id}_Unq.fasta"
    ${params.denoiser} -unoise3 "${sample_id}_Unq.fasta"  -zotus "${sample_id}_zotus.fasta" -tabbedout "${sample_id}_Unq_unoise3.txt" -minsize ${params.minAbundance}
    ${params.denoiser} -otutab ${upper_fasta} -zotus ${sample_id}_zotus.fasta -otutabout zotuTable.txt -mapout zmap.txt
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

  publishDir { params.illuminaDemultiplexed ? "07_blast" : "08_blast" }, mode: params.publishMode

  input:
    tuple val(sample_id), path(a), path(zotus_fasta), path(zotuTable), path(blast_db), path(custom_db)

  output:
    tuple val(sample_id), path("${sample_id}_blast_Result.tab")

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
    -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" \
    -perc_identity ${params.percentIdentity} -evalue ${params.evalue} \
    -best_hit_score_edge 0.05 -best_hit_overhang 0.25 \
    -qcov_hsp_perc ${params.qcov} -max_target_seqs ${params.maxQueryResults} \
    -query ${zotus_fasta} -out ${sample_id}_blast_Result.tab \
    -num_threads ${task.cpus}
  """
}

// make custom blast database for LULU curation
process lulu_blast {
  label 'blast'

  input:
    tuple val(sample_id), path(zotus_fasta)
  
  output:
    tuple val(sample_id), path('match_list.txt')

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
  label 'lulu'

  errorStrategy 'ignore'

  publishDir { params.illuminaDemultiplexed ? "08_lulu" : "09_lulu" }, mode: params.publishMode

  input:
    tuple val(sample_id), path(match_list), path(zotuTable)

  output:
    tuple path("curated_zotuTable.tab"), path("lulu_zotu_map.tab"), path("lulu_result_object.rds")

  script:
  """
  lulu.R ${params.luluMin}
  """
}

// retrieve one of the pre-trained insect models from https://github.com/shaunpwilkinson/insect#classifying-sequences
// because of previous sanity checks, we assume the value passed in `model` is a real one
process get_model {
  input:
    val(model)
  output:
    path('classifier.rds')

  exec:
    classifier = task.workDir / 'classifier.rds' 
    url = new URL(helper.insect_classifiers[model.toLowerCase()])
    url.withInputStream { stream -> classifier << stream }
}

// run insect classifier model
process insect_classify {
  label 'insect'

  publishDir { 
    p = params.assignTaxonomy ? "b" : ""
    return params.illuminaDemultiplexed ? "09${p}_insect_classified" : "10${p}_insect_classified" 
  }, mode: params.publishMode

  input:
    tuple path(zotus), path(classifier)

  output:
    tuple path('insect_classified.csv'), path('insect_settings.txt'), path('classifier.rds')

  script:
  """
  if [ "${classifier}" != "classifier.rds" ]; then
    mv ${classifier} classifier.rds
  fi
  insect.R ${zotus} classifier.rds \
    ${task.cpus} ${params.insectThreshold} \
    ${params.insectOffset} ${params.insectMinCount} \
    ${params.insectPing}
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

  // give example of what a demultiplexed FASTA file looks like
  if (params.demuxedExample) {
    helper.demuxed_example()
    exit(0)
  }

  // make sure the right version of single,paired,demultiplexed is passed
  if (!helper.file_exists(params.demuxedFasta) && params.single == params.paired) {
    if (!params.single) {
      println("One of either --single or --paired MUST be passed")
    } else {
      println("Only one of either --single or --paired may be passed")
    }
    exit(1)
  }

  // if denoiser is an executable, treat it as such
  // otherwise check to make sure it's a valid input
  if (helper.executable(params.denoiser)) {
    exec_denoiser = true
  } else {
    if (!(params.denoiser in ['usearch','usearch32','vsearch'])) {
      println("--denoiser must be either 'usearch', 'usearch32', 'vsearch', or a path to an executable (e.g., /opt/sw/bin/usearch64)")
      exit(1)
    }
  }

  // sanity check, blast database
  if (!helper.is_dir(params.blastDb)) {
    println("BLAST database must be specified either with the --blast-db argument")
    println("or using the \$BLASTDB environment variable. It must point to the directory")
    println("containing the `nt` database (do not include /nt in the path)")
    exit(1)
  }

  // make sure custom blast database is specified correctly
  if (params.customDbName != "NOTHING" || params.customDb != "NOTHING") {
    if (!helper.is_dir(params.customDb)) {
      println("Custom BLAST database must be specified as follows:")
      println("--custom-db-dir <path to custom BLAST db directory>")
      println("--custom-db <name of custom BLAST db (basename of .ndb, etc. files)>")
      println("example:")
      println("--custom-db-dir /storage/blast --custom-db test")
      exit(1)
    }
    if (params.customDbName == "NOTHING") {
      println("Name of custom BLAST db must be specified with the --custom-db option")
      exit(1)
    }
  }

  // another blast database sanity check
  if (helper.basename(params.blastDb) == helper.basename(params.customDb)) {
    println("Due to the vicissitudes of nextflow internality, the directory names")
    println("of the main and custom BLAST databases must be different.")
    println("As specified, both reside in directories called ${helper.basename(params.blastDb)}")
    exit(1)
  }

  // make sure insect parameter is valid: either a file or one of the pretrained models
  if (params.insect) {
    if (!helper.insect_classifiers[params.insect.toLowerCase()]) {
      if (!helper.file_exists(params.insect)) {
        println("Value passed to --insect must be one of the supported builtins or an RDS file")
        println("containing a trained insect classifier model.")
        println("See eDNAFlow.nf --help for supported builtin models")
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
include { taxonomy as assign_collapse_taxonomy  } from './modules/modules.nf'
include { taxonomy as assign_collapse_taxonomy_lulu  } from './modules/modules.nf'

workflow {
  // make sure our arguments are all in order
  check_params()

  vsearch = (params.vsearch || params.denoiser == 'vsearch')

  if (!helper.file_exists(params.demuxedFasta)) {
    if (params.single) {
      // here we load whatever was passed as the --reads option
      // if it's a glob, we get a list of files. if it's just one, we get just one
      // and we try to pull off something from the beginning to use as a sample ID. 
      Channel.fromPath(params.reads, checkIfExists: true) |
        map { [it.baseName.tokenize('_')[0],it] } |
        set { reads }
    } else {
      // if fwd and rev point to files that exists, just load them directly
      if ( helper.file_exists(params.fwd) && helper.file_exists(params.rev) ) {
        Channel.of(params.prefix) |
          combine(Channel.fromPath(params.fwd,checkIfExists: true)) | 
          combine(Channel.fromPath(params.rev,checkIfExists: true)) | 
          map { a,b,c -> [a,[b,c]] } |
          set { reads }
      } else {
        // otherwise make a glob pattern to find where the fwd/rev reads live
        // this may get a little complex, so there are two possible options:
        // files must have some way of stating which direction they are (e.e., R1/R2)
        // fwd/rev reads may optionally be in separate subdirectories but those must both be subdirectories
        // at the same level. So you can have /dir/R1 and /dir/R2, but you CAN'T have /dir1/R1 and /dir2/R2
        // as of now they must not be gzipped 
        if (params.fwd != "" && params.rev != "") {
          pattern = "${params.baseDir}/{${params.fwd},${params.rev}}/*{${params.r1},${params.r2}}*.fastq"
        } else {
          pattern = "${params.baseDir}/*{${params.r1},${params.r2}}*.fastq"
        }
  
        // load the reads and make sure they're in the right order because the
        // files will be sorted alphabetically. this is an extremely niche issue
        // because file are basically always R1/R2 but what if they're
        // forwards/backwards or something like that?

        // we also replace dashes with underscores in the sample ID because at least
        // vsearch will cut on that character and not treat it as a complete sample id
  
        // I'm not even certain the order matters, but I think it does because we 
        // send them to --file1 and --file2 of AdapterRemoval
        Channel.fromFilePairs(pattern) | 
          map { key,f -> [key.replaceAll("-","_"),f[0] =~ /${params.r1}/ ? [f[0],f[1]] : [f[1],f[0]]]  } | 
          map { key,f -> key == "" ? [params.prefix,f] : [key,f] } | 
          set { reads }
      }
    }

    // load barcodes
    barcodes = Channel.fromPath(params.barcode, checkIfExists: true)

    // if the sequences are already demultiplexed by illumina, we'll
    // process them separately, including optionally attempting to remove ambiguous tags
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

      // remove ambiguous tags, if specified
      if (params.removeAmbiguousTags) {
        reads_filtered_merged | 
          filter_ambiguous_tags | 
          set { reads_filtered_merged }
      } 

      // run the rest of the pipeline, including the primer mismatch check,
      // length filtering, and smashing together into one file
      reads_filtered_merged | 
        combine(barcodes) | 
        ngsfilter | 
        filter_length | 
        (vsearch ? relabel_vsearch : relabel_usearch) | 
        collectFile(name: "${params.prefix}_relabeled.fasta", storeDir: '05_relabeled_merged') |
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
            map { key, read1, read2 -> [key, [read1,read2]] } |
            set { reads }
        } else {
          // in single-end mode we can just split directly
          reads |
            splitFastq(by: params.splitBy, file: true) | 
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

      // run the rest of the pipeline, including demultiplexing, length filtering,
      // splitting, and recombination for dereplication
      reads_filtered_merged | 
        combine(barcodes) |
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
        // relabel to fasta
        (vsearch ? relabel_vsearch : relabel_usearch) | 
        // collect to single relabeled fasta
        collectFile(name: "${params.prefix}_relabeled.fasta", storeDir: '06_relabeled_merged') |
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
  Channel.of(params.prefix) | 
    combine(to_dereplicate) | 
    (vsearch ? derep_vsearch : derep_usearch) |
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
      map { sid, uniques, zotus, zotutable -> [sid,zotus] } | 
      lulu_blast | 
      combine(zotu_table) | 
      lulu
  }

  // run the insect classifier, if so desired
  // this should run in parallel with the blast & lulu processes
  if (params.insect) {
    // dereplicate returns a tuple, but we only need the zotus fasta
    dereplicated | 
      map { sid, uniques, zotus, zotutable -> zotus } | 
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
    zotus | 
      combine(classifier) | 
      insect_classify
  }


  // run taxonomy assignment/collapse script if so requested
  if (params.assignTaxonomy && !params.skipBlast) {
    // here we grab the blast result
    blast.out | 
      map { sid, blast_result -> blast_result } | 
      set { blast_result }


    // then we smash it together with the blast results 
    // and run the taxonomy assignment/collapser script
    zotu_table |
      combine(blast_result) |
      combine(Channel.of('uncurated')) | 
      assign_collapse_taxonomy

    if (!params.skipLulu) {
      // run the taxonomy assignment for lulu-curated zotus
      lulu.out | 
        map { zotutable, zotu_map, result_object -> zotutable } | 
        combine(blast_result) |
        combine(Channel.of('curated')) | 
        assign_collapse_taxonomy_lulu
    }
  }
}
