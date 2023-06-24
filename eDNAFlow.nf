#!/usr/bin/env nextflow20
nextflow.enable.dsl=2

/* some global variables */
exec_denoiser = false

// TODO: dont' use skipDemux. instead skip demultiplexing just if demuxedFasta is an existing file
// TODO: incorporate this method of reading single or paired reads: https://github.com/nextflow-io/nextflow/issues/236#issuecomment-314018546
// TODO: do insect assignment. here's a docker image: docker://mahsamousavi/insect:2019
// TODO: make the taxonomy collapser script more legible and easier to use, have it cache taxdump, etc.
// TODO: since we stopped doing the illumina_sample annotation, we probably still need to do the
//       ambiguous tag removal step somewhere
// TODO: consider putting back in the stats (sequence counts, sample names, etc.)
// TODO: rather than putting parameters in filenames, store parameter values in a text file in each output dir

// gets an environment variable and returns blank instead of null if it's not set

// FastQC to check the initial quality of raw reads
process initial_fastqc {
  label 'fastqc'

  publishDir "00_fastqc_initial", mode: params.publishMode

  input:
    tuple val(sample_id), path(read) 

  output:
    path('*_fastqc.{zip,html}')

  script:
  if( read instanceof Path ) {
    // single end
    """
    fastqc -q ${read}
    """
  } else {
    // paired end
    """
    fastqc -q ${read[0]} ${read[1]}
    """
  }
}

// multiqc combines the output of multiple fastqc reports
// we do this for the non-demultiplexed samples
process initial_multiqc {
  label 'multiqc'

  publishDir '00_fastqc_initial', mode: params.publishMode

  input:
    path(fastqc_files)
  output: 
    tuple path('multiqc_report.html'), path('multiqc_data/*')

  script:
  """
  multiqc .
  """
}

// FastQC to check quality after filtering and merging
process merged_fastqc {
  label 'fastqc'

  publishDir "00_fastqc_filtered_merged", mode: params.publishMode

  input:
    tuple val(sample_id), path(read) 

  output:
    path('*_fastqc.{zip,html}') 

  script:
  if( read instanceof Path ) {
    // single end
    """
    fastqc -q ${read}
    """
  } else {
    // paired end
    """
    fastqc -q ${read[0]} ${read[1]}
    """
  }
}

process merged_multiqc {
  label 'multiqc'

  publishDir '00_fastqc_filtered_merged', mode: params.publishMode

  input:
    path(fastqc_files)
  output: 
    tuple path('multiqc_report.html'), path('multiqc_data/*')

  script:
  """
  multiqc .
  """
}


// trim and (where relevant) merge paired-end reads
// this gets called if the files were already demultiplexed
// TODO: I'm not sure we need two entirely separate processes for this
process filter_merge_demultiplexed  {
  label 'adapterRemoval'

  publishDir '01_demultiplexed_filter_merge', mode: params.publishMode

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

// annotate sequences that have already been demultiplexed by illumina
// we probably don't actually need this but it's being left in for now
process annotate_demultiplexed {
  label 'obitools'

  /* TODO: allow customization of how the sample IDs are interpreted */
  /* for now we will use the sample id returned by fromFilePairs */

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

// filter, trim, and (where relevant) merge reads
// this gets called for raw reads that were *not* demultiplexed by the sequencer
process filter_merge {
  label 'adapterRemoval'

  cpus { params.maxCpus }

  publishDir "01_filter_merge", mode: params.publishMode

  input:
  tuple val(sample_id), path(read) 

  output:
  tuple val(sample_id), path('*_QF.fastq') 

  script:

  if( read instanceof Path ) {   
    // single end
    """
    AdapterRemoval --threads ${task.cpus} --file1 ${read} \
      --trimns --trimqualities \
      --minquality ${params.minQuality} \
      --basename ${sample_id}

    mv ${sample_id}.truncated ${sample_id}_QF.fastq
    """
  } else {  
    // if reads are paired-end then merge 
    """
    AdapterRemoval --threads ${task.cpus} --file1 ${read[0]} --file2 ${read[1]} \
      --collapse --trimns --trimqualities \
      --minquality $params.minQuality \
      --minalignmentlength ${params.minAlignLen} \
      --basename ${sample_id}

    mv ${sample_id}.collapsed ${sample_id}_QF.fastq  
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
  ngsfilter --uppercase -t ${barcode} -e ${params.primerMismatch} -u "orphan.fastq" ${read} > "${sample_id}_${barcode.baseName}_QF_Dmux.fastq"
  """
}


// combine outputs from (possible) multiple barcode files, filter by length
process filter_length {
  label 'obitools'

  publishDir "03_length_filtered"

  input: 
    tuple val(sample_id), val(barcode_files), path(fastq_files) 
  
  output:
    tuple val(sample_id), path('*_QF_Dmux_minLF.fastq') 

  script:
  // if we're already demultiplexed we probably don't have the forward_tag and reverse_tag annotations
  def p = params.illuminaDemultiplexed ? "" : "-p 'forward_tag is not None and reverse_tag is not None'" 
  """
  cat ${fastq_files} > "${sample_id}_QF_Dmux.fastq" 
  obigrep --uppercase -l ${params.minLen} ${p} "${sample_id}_QF_Dmux.fastq" > "${sample_id}_QF_Dmux_minLF.fastq"
  """
}



// for non-demultiplexed runs, split the annotated reads file by samples
process split_samples {
  label 'obitools'

  publishDir "04_splitSamples_${sample_id}", mode: params.publishMode

  input:
    tuple val(sample_id), path(fastqs) 
  
  output:
    /*tuple val(sample_id),*/ path "$sample_id/*.fastq"

  script:
  // TODO: we probably don't need this anymore since we're not doing this annotation now
  def sample_tag = params.illuminaDemultiplexed ? "illumina_sample" : "sample"
  """
  mkdir ${sample_id}
  obisplit --uppercase -t ${sample_tag} -u "noSampleID.fastq" $fastqs
  mv *.fastq ${sample_id}
  mv ${sample_id}/$fastqs ..
  mv ${sample_id}/noSampleID.fastq  noSampleID.fastq.ignore
  """
}

// relabel files for vsearch
process relabel_vsearch {
  label 'vsearch'

  cpus { params.maxCpus }

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

    # for files in \${fastqs}; do
    #   label=\$(echo \$files | cut -d '/' -f 3 | cut -d '.' -f 1)
    #   usearch -fastx_relabel \$files -prefix \${label}. -fastqout \${label}.relabeled.fastq 
    # done

    # for files in *.relabeled.fastq; do
    #   name=\$(echo \$files | cut -d '/' -f '2' | cut -d '.' -f 1)
    #   echo \${name} >> CountOfSeq.txt
    #   grep "^@\${name}" \$files | wc -l >> CountOfSeq.txt
    # done 

    # cat *.relabeled.fastq > "${sample_id}_QF_Dmux_minLF_relabeled4Usearch.fastq"

    # usearch -fastx_get_sample_names *_relabeled4Usearch.fastq -output sample.txt

    # usearch -fastq_filter *_relabeled4Usearch.fastq -fastaout ${sample_id}_upper.fasta


    """
  } else if (exec_denoiser) {
    // TODO: will this work? recall that it will be executed in the context of the usearch docker image
    //       this is why calling vsearch like this didn't work
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

  // give us all the processors, since we should be the only ones on this
  cpus { params.maxCpus }

  publishDir { params.illuminaDemultiplexed ? "05_derep_vsearch" : "06_derep_vsearch" }, mode: params.publishMode

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

  publishDir { params.illuminaDemultiplexed ? "05_derep_usearch" : "06_derep_usearch" }, mode: params.publishMode

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

  // get all the cpus and memory we can
  cpus { params.maxCpus }
  memory { params.maxMemory }

  publishDir { params.illuminaDemultiplexed ? "07_blast" : "08_blast" }, mode: params.publishMode

  input:
    tuple val(sample_id), path(a), path(zotus_fasta), path(zotuTable), path(blast_db), path(custom_db)

  output:
    tuple val(sample_id), path("${sample_id}_blast_Result.tab"), path(zotuTable), path("match_list.txt") 

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

  publishDir { params.illuminaDemultiplexed ? "08_lulu" : "09_lulu" }, mode: params.publishMode

  input:
    tuple val(sample_id), path(a), path(zotuTable), path(match_list) 

  output:
    tuple path("curated_zotuTable.tab"), path("lulu_zotu_map.tab"), path("lulu_result_object.rds")

  script:
  """
  lulu.R ${params.luluMin}
  """
}


// taxonomy assignment/collapse to least common ancestor (LCA)
process assign_collapse_taxonomy {
  label 'lca_python3'
    publishDir { params.illuminaDemultiplexed ? '09_collapsed_taxonomy' : '10_collapsed_taxonomy' }, mode: params.publishMode

    input:
      tuple path(zotu_table), path(lulu_zotu_table), path(blast_result) 

    output:
      tuple path("intermediate_table.tab"), path("taxonomy_collapsed.tab"), path("intermediate_table_lulu.tab"), path("taxonomy_collapsed_lulu.tab") // path("${params.lcaOutput}_qCov${params.lcaQcov}_id${params.lcaPid}_diff${params.lcaDiff}.tab") 


    script:
    """
    runAssign_collapsedTaxonomy.py ${zotu_table} ${blast_result} ${params.lcaQcov} ${params.lcaPid} ${params.lcaDiff} taxonomy_collapsed.tab
    mv interMediate_res.tab intermediate_table.tab
    runAssign_collapsedTaxonomy.py ${lulu_zotu_table} ${blast_result} ${params.lcaQcov} ${params.lcaPid} ${params.lcaDiff} taxonomy_collapsed_lulu.tab
    mv interMediate_res.tab intermediate_table_lulu.tab
    """
}  

def basename(f) {
  (new File(f)).baseName
}

def file_exists(f) {
  (new File(f)).exists()
}

def is_dir(d) {
  (new File(d)).isDirectory()
}

def executable(f) {
  (new File(f)).canExecute()
}

def get_env(n) {
  v = System.getenv(n)
  return v ? v : ""
}

def usage() {
  log.info """
  Usage: eDNAFlow.nf [options]

  General options:
    --prefix [str]             Project prefix, applied to output filenames as "sample ID" 
                               (default: 'seq')
    --barcode [file]           (required) Barcode file. Format must match OBITools requirements
                               (see https://pythonhosted.org/OBITools/scripts/ngsfilter.html)
                               To denote multiple barcode files, this may be a glob, but it must
                               be surrounded by quotations (e.g. 'barcode*.tab'). 
                               For previously-demultiplexed datasets, use ':' for barcode pairs
                               Primer sequences are still used for primer-mismatch comparisons
    --publish-mode [str]       Specify how nextflow places files in output directories 
                               (default: symlink)
    --fastqc                   Output FastQC reports for pre and post filter/merge steps 

  For single-end sequencing runs:
    --single                   Specify single-ended sequencing run (required)
    --reads [file/glob]        Location of sequencing read(s). For runs that have already been 
                               demultiplexed by the sequencer, this may be a glob (e.g., '*.fastq')

  For paired-end sequencing runs:  
    --paired                   Specify paired-end sequencing run (required)

  To specify the location of paired-end sequence files, the following file patterns are followed:
    <base-dir>/*<r1>|<r2>*.fastq
    <base-dir>/<fwd>|<rev>/*<r1>|<r2>*.fastq

  Options:
    --base-dir [dir]           Base directory of sequence read subdirectories (default: .) 
    --fwd [dir]                (Optional) forward reads directory (default: ${params.fwd})
    --rev [dir]                (Optional) reverse reads directory (default: ${params.rev})
    --r1 [str]                 Pattern distinguishing forward read files (default: ${params.r1})
    --r2 [str]                 Pattern distinguishing reverse read files (default: ${params.r2})

  Length, quality, and merge settings:
    --min-quality [num]        Minimum Phred score for sequence retention 
                               (default: ${params.minQuality})
    --min-align-len [num]      Minimum sequence overlap when merging forward/reverse reads
                               (default: ${params.minAlignLen}
    --min-len [num]            Minimum overall sequence length (default: ${params.minLen})

  BLAST settings (one or more of --blast-db or --custom-db is required):
    --blast-db [dir]           Location of local BLAST nucleotide (nt) database *directory*
                               (do NOT include the final "/nt")
                               May be specified using the \$BLASTDB environment variable
                               (current value: ${get_env("BLASTDB")})
    --custom-db [dir]          (optional) Path to custom BLAST database *directory*
                               Passing --custom-db and --blast-db values that point to directories
                               with the same name (e.g., /dir1/blast and /dir2/blast) will result 
                               in an error!
    --custom-db-name [str]     Name of custom BLAST database (i.e., basename of .ndb, etc. files)
    --blast-task [str]         Set blast+ task (default: "blastn")
    --max-query-results [num]  Maxmimum number of BLAST results to return per zOTU (default: 10)
    --percent-identity [num]   Minimum percent identity of matches to report (default: 95)
    --evalue [num]             Expectation value threshold for saving hits (default: 0.001)
    --qcov [num]               Percent query coverage per hsp (default: 100)

  Taxonomy assignment:
    --assign-taxonomy          Perform final taxonomy assignment & LCA collapse 
    --lca-qcov [num]           Minimum query coverage for LCA taxonomy assignment (default: 100)
    --lca-pid [num]            Minimum percent identity for LCA taxonomy assignment (default: 97)
    --lca-diff [num]           The difference between percent identities (when query coverage is
                               identical) where species-level taxonomy is retained (default: 1)
    --insect [str]             Perform taxonomy assignment using insect and specify which classifier
                               model to use (not yet supported)

  Demultiplexing and sequence matching:
    --illumina-demultiplexed   Sequencing run has already been demultiplexed
    --remove-ambiguous-tags    Removes reads with ambiguous tags in the header (i.e., not A,G,C,T)
                               (only applies to previously-demultiplexed runs with tags in header)
    --demuxed-fasta [file]     Skip demultiplexing step and use supplied FASTA 
                               (must be in usearch/vsearch format)
    --demuxed-example          Spit out example usearch/vsearch demultiplexed FASTA format
    --demux-only               Stop after demultiplexing and splitting raw reads
    --primer-mismatch          Allowed number of mismatched primer bases 
                               (default: ${params.primerMismatch})

  Denoising and zOTU inference:
    --min-abundance [num]      Minimum zOTU abundance; zOTUs below threshold will be discarded
                               (default: ${params.minAbundance}) 
    --denoiser [tool/path]     Sets the tool used for denoising & chimera removal
                               accepted options: usearch/usearch32 (equivalent), vsearch, 
                               path to usearch64 executable (default: usearch)

  LULU zOTU curation:
    --lulu-min [num]           Minimum threshold of sequence similarity to consider zOTUs as spurious.
                               Choose higher values when using markers with lower genetic variation 
                               and/or few expected PCR and sequencing errors. (default: ${params.luluMin})

  Resource allocation:
    --max-memory [str]         Maximum memory available to all nextflow processes, e.g., '8.GB' (default: ${params.maxMemory})
    --max-cpus [num]           Maximum cores available to all nextflow processes default: ${params.maxCpus})
    --max-time [str]           Maximum time allocated to each pipeline process, e.g., '2.h' (default: ${params.maxTime})

  Singularity options:
    --bind-dir [str]           Space-separated list of directories to bind within singularity images
                               (must be surrounded by quotations if more than one directory)
                               Note: singularity will attempt to auto-bind all provided host paths
                               so this option may not be necessary, but try it if you're getting 
                               weird "file not found" types of errors
    --singularity-cache [str]  Location to store singularity images. May also be specified
                               with the environment variable \$NXF_SINGULARITY_CACHEDIR.
                               (current value: ${get_env("NXF_SINGULARITY_CACHEDIR")})
  """.stripIndent()
}

def demuxed_example() {
  log.info """
  >sample1.1
  AGCGTCCGATGACTGACTGACTAGCT
  >sample1.2
  TACGTACGATCGACGAGTCTACGACTACTGAC
  >sample1.3
  TGACTGATCGTACTATCAGAGCTATCATCGACTATCATCGATC
  >sample2.1
  ATCGTACTACTAGCGACGAGTCATCACGACGTACTAGTCGA
  >sample2.2
  CATGCGACGTACGTACTATCATCATCGAGCAGCTATATATCGATGGTACTAGCTGAC
  >sample2.3
  TGACTGATCGTACTATCAGAGCTATCATCGACTATCATCGATC
  >sample3.1
  AGCGTCCGATGACTGACTGACTAGCT
  >sample3.2
  ATCGTACTACTAGCGACGAGTCATCACGACGTACTAGTCGA
  >sample3.3
  CATGCGACGTACGTACTATCATCATCGAGCAGCTATATATCGATGGTACTAGCTGAC
  """.stripIndent()
}

def check_params() {
// show help message and bail
  if (params.help) {
    usage()
    if (params.debug) {
      println("\n\n\n")
      println(params)
    }
    exit(0)
  }

  if (params.demuxedExample) {
    demuxed_example()
    exit(0)
  }

  if (params.single == params.paired) {
    if (!params.single) {
      println("One of either --single or --paired MUST be passed")
    } else {
      println("Only one of either --single or --paired may be passed")
    }
    exit(1)
  }

  if (executable(params.denoiser)) {
    exec_denoiser = true
  } else {
    if (!(params.denoiser in ['usearch','usearch32','vsearch'])) {
      println("--denoiser must be either 'usearch', 'usearch32', 'vsearch', or a path to an executable (e.g., /opt/sw/bin/usearch64)")
      exit(1)
    }
  }

  if (!is_dir(params.blastDb)) {
    println("BLAST database must be specified either with the --blast-db argument")
    println("or using the \$BLASTDB environment variable. It must point to the directory")
    println("containing the `nt` database (do not include /nt in the path)")
    exit(1)
  }

  if (params.customDbName != "NOTHING" || params.customDb != "NOTHING") {
    if (!is_dir(params.customDb)) {
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

  if (basename(params.blastDb) == basename(params.customDb)) {
    println("Due to the vicissitudes of nextflow internality, the directory names")
    println("of the main and custom BLAST databases must be different.")
    println("As specified, both reside in directories called ${basename(params.blastDb)}")
    exit(1)
  }
}


workflow {

  // make sure our arguments are all in order
  check_params()

  if (!params.skipDemux) {
    if (params.single) {
      // here we load whatever was passed as the --reads option
      // if it's a glob, we get a list of files. if it's just one, we get just one
      // and we try to pull off something from the beginning to use as a sample ID. 
      // TODO: customize how to specify sample IDs
      Channel.fromPath(params.reads, checkIfExists: true) |
        map { it -> [it.baseName.tokenize('_')[0],it] } |
        set { reads }
    } else {
      // if fwd and rev point to files that exists, just load them directly
      if ( file_exists(params.fwd) && file_exists(params.rev) ) {
        Channel.of(params.prefix) |
          combine(Channel.fromPath(params.fwd,checkIfExists: true)) | 
          combine(Channel.fromPath(params.rev,checkIfExists: true)) | 
          map { a,b,c -> [a,[b,c]] } |
          set { reads }
      }
      else {
        // otherwise make a glob pattern to find where the fwd/rev reads live
        // this may get a little complex, so there are two possible options:
        // files must have some way of stating which direction they are (e.e., R1/R2)
        // fwd/rev reads may optionally be in separate subdirectories but those must both be subdirectories
        // at the same level. So you can have /dir/R1 and /dir/R2, but you CAN'T have /dir1/R1 and /dir2/R2
        // as of now they must not be gzipped (though TODO: we can do the gunzipping, that's easy)
        if (params.fwd != "" && params.rev != "") {
          pattern = "${params.baseDir}/{${params.fwd},${params.rev}}/*{${params.r1},${params.r2}}*.fastq"
        } else {
          pattern = "${params.baseDir}/*{${params.r1},${params.r2}}*.fastq"
        }
  
        // load the reads and make sure they're in the right order because the
        // files will be sorted alphabetically this is an extremely niche issue
        // because file are basically always R1/R2 but what if they're
        // forwards/backwards or something like that?
  
        // I'm not even certain the order matters, but I think it does because we 
        // send them to --file1 and --file2 of AdapterRemoval
        // TODO: figure out how to customize sample IDs here too, maybe
        Channel.fromFilePairs(pattern) | 
          map { key,f -> [key,f[0] =~ /${params.r1}/ ? [f[0],f[1]] : [f[1],f[0]]]  } | 
          set { reads }
      }
    }

    // load barcodes
    barcodes = Channel.fromPath(params.barcode, checkIfExists: true)

    // if the sequences are already demultiplexed by illumina, we'll
    // annotate them with their sample names and optionally filter out
    // sequences with ambiguous tags. then we recombine them into one forward
    // and one reverse file and continue with those as our reads
      if (params.illuminaDemultiplexed) {
        // run the first part of the pipeline for sequences that have already
        // been demultiplexed by the sequencer

        if (params.fastqc) {
          reads | 
            initial_fastqc | 
            flatten | 
            collect | 
            initial_multiqc
        }

        reads |
          filter_merge_demultiplexed |
          set { reads_filtered_merged }

        if (params.fastqc) {
          reads_filtered_merged | 
            merged_fastqc | 
            flatten | 
            collect | 
            merged_multiqc
        }

        reads_filtered_merged | 
          /* annotate_demultiplexed |  */
          combine(barcodes) | 
          ngsfilter | 
          groupTuple | 
          filter_length | 
          (params.denoiser == 'vsearch' ? relabel_vsearch : relabel_usearch) | 
          collectFile(name: "${params.prefix}_relabeled.fasta", storeDir: '06_relabeled_merged') |
          set { to_dereplicate } 
      } else {
        // do initial fastqc step
        if (params.fastqc) {
          initial_fastqc(reads)
        }

        // do quality filtering and/or paired-end merge
        reads | 
          filter_merge | 
          set { reads_filtered_merged }

        // do after-filtering fastqc step
        if (params.fastqc) {
          merged_fastqc(reads_filtered_merged)
        }

        reads_filtered_merged | 
          combine(barcodes) |
          ngsfilter | 
          groupTuple | 
          filter_length |
          split_samples | 
          flatten | 
          map { it -> [it.baseName, it] } | 
          (params.denoiser == 'vsearch' ? relabel_vsearch : relabel_usearch) | 
          collectFile(name: "${params.prefix}_relabeled.fasta", storeDir: '07_relabeled_merged') |
          set { to_dereplicate }
      }
  } else {
    // here we've already demultiplexed and relabeled sequences 
    // (presumably from an earlier run of the pipeline), so we can jump to here

    Channel.of(params.demuxedFasta) | 
      combine(Channel.fromPath(params.demuxedFasta)) |
      set { to_dereplicate }
  }

  // build the channel, run dereplication, and set to a channel we can use again
  Channel.of(params.prefix) | 
    combine(to_dereplicate) | 
    (params.denoiser == 'vsearch' ? derep_vsearch : derep_usearch) |
    set { dereplicated }

  // make blastdb and customDbName a path channel so singularity will automount it properly
  blast_db = Channel.fromPath(params.blastDb, type: 'dir')
  custom_db = Channel.fromPath(params.customDb, type: 'dir', checkIfExists: params.customDb != "NOTHING")

  // run blast query and lulu curation
  dereplicated | 
    combine(blast_db) | 
    combine(custom_db) | 
    blast | 
    lulu


  // run taxonomy assignment/collapse script if so requested
  if (params.assignTaxonomy) {
    // here we grab the zotu table from our dereplication step
    // TODO: also include lulu-curated output
    dereplicated | 
      map { it -> it[3] } | 
      set { zotu_table }

    // here we grab the blast result
    blast.out | 
      map { it -> it[1] } | 
      set { blast_result }

    lulu.out | 
      map { it -> it[0] } | 
      set { lulu_zotu_table }

    // then we smash them together and run the taxonomy assignment/collapser script
    zotu_table |
      combine(lulu_zotu_table) | 
      combine(blast_result) |
      assign_collapse_taxonomy
  }

}
