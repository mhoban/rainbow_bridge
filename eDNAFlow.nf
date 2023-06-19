#!/usr/bin/env nextflow20
nextflow.enable.dsl=2

def usage() {
  println("Usage:")
  println("")
  println("Examples of basic commands to run the pipeline on local machine:")
  println("")
  println("For single-end run:")
  println("nextflow run eDNAFlow.nf --reads 'file.fastq' --barcode 'bc_*.txt' --blast_db 'path2/LocalGenbankDatabase/nt' [OPTIONS] ")
  println("")
  println("For paired-end run:")
  println("nextflow run eDNAFlow.nf --barcode 'pe_bc*' --blast_db 'Path2TestBlastDataset/file.fasta' --custom_db 'path2/customDatabase/myDb' [OPTIONS]")
  println("")
  println("For running LCA taxonomy assignment script")
  println("nextflow run eDNAFlow.nf --taxonomyAssignment --zotuTable \"path2/curatedOruncurated_ZotuTable_file\" --blastFile \"path2/blastResult_file\" --lca_output \"my_lca_result\" [OPTIONS]")
  println("")
  println("Mandatory arguments if sequences are NOT demultiplexed")
  println("  --reads [file]                 Input data name (must be surrounded with quotes);")
  println("                                 do NOT specify this option if reads are paired-end as they get identified automatically by default;")
  println("                                 reads must be unzipped. The paired-end file name must end with _R1.fastq & _R2.fastq")
  println("  --barcode [file]               Barcode file name; barcode file format must match OBITools requirement;")
  println("                                 if multiple barcode files exist (e.g. bc_1.txt; bc_2.txt) it can be specified like this: 'bc*.txt'  ")
  println("  ")
  println("At least one of the below databases must be specified.")
  println("  --blast_db [dir]               Absolute path to local nt databse ")
  println("  --custom_db [dir]              Absolute path to custom database")
  println("")
  println("Mandatory arguments if sequences are demultiplexed")
  println("  --skipDemux [bool]")
  println("  --demuxedInput [file]          A fasta file holding all the demultiplexed sequences;")
  println("                                 Format of the sample identifier must match USEARCH requirements")
  println("")
  println("At least one of the below databases must be specified.")
  println("  --blast_db [dir]               Absolute path to local nt databse ")
  println("  --custom_db [dir]              Absolute path to custom database")
  println("")
  println("Mandatory and optional arguments for running LCA taxonomy assignment script")
  println("  --taxonomyAssignment [bool]")
  println("  --zotuTable [file]             A raw or LULU curated Zotu, OTU, or ASV table file; file must exist in the same directory as eDNAFlow;")
  println("  --blastFile [file]             Blast result file; file must exist in the same directory as eDNAFlow")
  println("                                 For file format requirements of --zotuTable & --blastFile check out eDNAFlow GitHub page")
  println("")
  println("Optional")
  println("  --lca_qcov    [num]            Percent of query coverage; Default is ${params.lca_qcov}")
  println("  --lca_pid     [num]            Percent of identity; Default is ${params.lca_pid}")
  println("  --lca_diff    [num]            The difference (Diff) between % identities of two hits when their qCov is equalDiff; Default is ${params.lca_diff}")
  println("  --lca_output  [string]         Output file name; Default is ${params.lca_output}")
  println("")
  println("Skipping                         Skip any of the mentioned steps")
  println("  --skipDemux [bool]             If this option is set, then --demuxedInput [file] must be provided")
  println("  --skipFastqc [bool]")
  println("")
  println("Parameters to run eDNAFlow on Cloud/HPC")
  println("  -profile [string]              Currently can choose between \"nimbus\" (can be used if user has access to more memory i.e. cloud or HPC),")
  println("                                 \"zeus\" and \"magnus\" (it's specific to users who have access to ZEUS/Magnus - high-throughput HPC clusters at the Pawsey Supercomputing Centre)  ")
  println("                                 e.g. -profile nimbus")
  println("  ")
  println("  --bindDir [dir]                If you run eDNAFlow on Cloud or HPC, you will need to specify this option")
  println("                                 On HPC, it usually will be /scratch and/or /group. On Cloud, it could be your mounted volume.")
  println("                                 e.g. --bindDir \"/scratch\"")
  println("                                 If you need to mount more than one directory put space in between e.g. --bindDir \"/scratch /group\"")
  println("")
  println("General optional parameters")
  println("  --help                         Show this help message")
  println("  --publish_dir_mode [string]    Choose between symlink , copy, link. Default is ${params.publish_dir_mode}")
  println("  --singularityDir [dir]         Directory where singularity images will be stored")
  println("")
  println("Demultiplexing")
  println("  --onlyDemux [bool]")
  println("  --primer_mismatch [num]        Number of mismatch allowed for matching primers; Default is ${params.primer_mismatch}")
  println("")
  println("Quality filtering / Merging")
  println("  --minQuality [num]             The minimum Phred quality score to apply for quality control of raw sequences; Default is ${params.minQuality}")
  println("  --minAlignLeng [num]           The minimum alignment length for merging read1 and read2; Default is ${params.minAlignLeng}")
  println("  --minLen [num]                 The minimum length allowed for sequences; Default is ${params.minLen}")
  println("")
  println("ZOTU formation")
  println("  --minsize [num]                The minimum abundance; input sequences with lower abundances are removed; Default is ${params.minsize}")
  println("")
  println("Blast parameters")
  println("  --blast_task [string]	        Blast task to be performed; Default is '${params.blast_task}'; but can be set to 'megablast' if required   ")
  println("  --maxTarSeq [num]              The maximum number of target sequences for hits per query to be returned by Blast; Default is ${params.maxTarSeq}")
  println("  --perc_identity [num]          Percentage of identical matches; Default is '${params.perc_identity}'")
  println("  --evalue [num]                 Expected value for saving blast hits; Default is '${params.evalue}'")
  println("  --qcov [num]                   The percent of the query that has to form an alignment against the reference to be retained;")
  println("                                 Higher values prevent alignments of only a short portion of the query to a reference; Default is '${params.qcov}'")
  println("  ")
  println("Choice of USEARCH32 vs USEARCH64 ")
  println("  --mode [str]                   Default is '${params.mode}'; for running with 64 version the mode has to be set to --mode 'usearch64'")
  println("                                 and below --search64 option has to be specified as well; can set to --mode 'vsearch'  ")
  println("  --usearch64 [dir]              Full path to where usearch64 bit version is stored locally")
  println("")
  println("LULU")
  println("  --lulu [file]                  An R script to run post-clustering curation with default settings of LULU;")
  println("                                 This file has been provided and must be present in the same directory as other scripts;")
  println("                                 by default eDNAFlow will be looking for this file in the same directory where eDNAFlow.nf is. ")
  println("  --minMatch_lulu [num]          Default is '${params.minMatch_lulu}'; A minimum threshold (minimum_match) of sequence similarity")
  println("                                 for considering any OTU as an error of another. ")
  println("                                 This setting should be adjusted so higher threshold is employed for genetic markers with little variation")
  println("                                 and/or few expected PCR and sequencing errors (See LULU paper). ")
  println("")
  println("Other options")
  println("  --max_memory [str]             Memory limit for each step of pipeline. e.g. --max_memory '8.GB'. Default: '${params.max_memory}'")
  println("  --max_time [str]               Time limit for each step of the pipeline. e.g. --max_time '2.h'. Default: '${params.max_time}'")
  println("  --max_cpus [str]               Maximum number of CPUs to use for each step of the pipeline. e.g. --max_cpus '1' Default: '${params.max_cpus}'")
}

// FastQC to check the initial quality of raw reads
process initial_fastqc {
  label 'fastqc'

  publishDir "00_fastQC_${qc_step}_${sample_id}", mode: params.publish_dir_mode

  input:
  tuple val(sample_id), path(read) 
  val qc_step

  output:
  tuple val(sample_id), path('*_fastqc.{zip,html}')

  when:
  !params.skipFastqc

  script:

  if( read instanceof Path ) {
    """
      fastqc -q ${read}
    """
  } else {
    """
      fastqc -q ${read[0]} ${read[1]}
    """
  }
}

// FastQC to check quality after filtering and merging
process merged_fastqc {
  label 'fastqc'

  publishDir "00_fastQC_${qc_step}_${sample_id}", mode: params.publish_dir_mode

  input:
  tuple val(sample_id), path(read) 
  val qc_step

  output:
  tuple val(sample_id), path('*_fastqc.{zip,html}') 

  when:
  !params.skipFastqc

  script:

  if( read instanceof Path ) {
    """
      fastqc -q ${read}
    """
  } else {
    """
      fastqc -q ${read[0]} ${read[1]}
    """
  }
}

 
// annotate sequences that have already been demultiplexed by illumina
process annotate_demultiplexed {
  label 'obitools'

  /* TODO: allow customization of how the sample IDs are interpreted */
  /* TODO: annotated the fwd/rev pair in parallel (can't do using the current docker img) */
  /* for now we will use the sample id returned by fromFilePairs */

  publishDir '00_annotated', mode: params.publish_dir_mode

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("*${params.r1}*_annotated.fastq"), path("*${params.r2}*_annotated.fastq")
    /* path "*${params.r1}*_annotated.fastq", emit: fwd */
    /* path "*${params.r2}*_annotated.fastq", emit: rev */

  script:
    def (f1,f2) = reads
    """
    # sample1=\$(basename ${f1} | cut -f1 -d_)
    # sample2=\$(basename ${f2} | cut -f1 -d_)
    # optionally remove ambiguous tags
    if ${params.remove_ambiguous_tags}; then
      # annotate sample names
      obiannotate --uppercase -S illumina_sample:"'${sample_id}'" ${f1} | obigrep --uppercase -D ':[ACGT]+\\+[ACGT]+\$' > "\$(basename ${f2} .fastq)_annotated.fastq"
      obiannotate --uppercase -S illumina_sample:"'${sample_id}'" ${f2} | obigrep --uppercase -D ':[ACGT]+\\+[ACGT]+\$' > "\$(basename ${f2} .fastq)_annotated.fastq"
    else
      # annotate sample names
      obiannotate --uppercase -S illumina_sample:"'${sample_id}'" ${f1} > "\$(basename ${f1} .fastq)_annotated.fastq"
      obiannotate --uppercase -S illumina_sample:"'${sample_id}'" ${f2} > "\$(basename ${f2} .fastq)_annotated.fastq"
    fi
    """
}

// Process 1: Quality filtering of raw reads
process filter_merge {
  label 'adapterRemoval'

  cpus { params.max_cpus / 2 }

  publishDir "01_a_quality_Filtering_${sample_id}", mode: params.publish_dir_mode

  input:
  tuple val(sample_id), path(read) 

  output:
  tuple val(sample_id), path('*_QF.fastq') 

  script:

  if( read instanceof Path ) {   
    """
      AdapterRemoval --threads ${task.cpus} --file1 ${read} \
      --trimns --trimqualities \
      --minquality ${params.minQuality} \
      --basename ${sample_id}

    mv ${sample_id}.truncated ${sample_id}_QF.fastq

      """
  }

// if reads are paired-end then merge 
  else {  
    """
      AdapterRemoval --threads ${task.cpus} --file1 ${read[0]} --file2 ${read[1]} \
      --collapse --trimns --trimqualities \
      --minquality $params.minQuality \
      --minalignmentlength ${params.minAlignLeng} \
      --basename ${sample_id}

    mv ${sample_id}.collapsed ${sample_id}_QF.fastq  
      """
  }
}

// Process 2: Assigning reads to samples/demultiplexing
process demultiplex {
  label 'obitools'

  publishDir "02_assigned_dmux_${sample_id}_${barcode.baseName}", mode: params.publish_dir_mode

  input:
    tuple val(sample_id), path(read), path(barcode) 


  output:
    tuple val(sample_id), val("${barcode.baseName}"), path("*_Dmux.fastq") 

  script:
  """
  ngsfilter --uppercase -t ${barcode} -e ${params.primer_mismatch} -u "orphan.fastq" ${read} > "${sample_id}_${barcode.baseName}_QF_Dmux.fastq"
  """
}


// Process 3: Cat multiple file after demux, Filtering sequences shorter than min length and single barcoded reads
process filter_length {
  label 'obitools'

  publishDir "03_Length_filtered_${sample_id}", mode: params.publish_dir_mode

  input: 
    tuple val(sample_id), val(barcode_files), path(fastq_files) 
  
  output:
    tuple val(sample_id), path('*_QF_Dmux_minLF.fastq') //into lenFilt_ch


  script:
  """
  cat ${fastq_files} > "${sample_id}_QF_Dmux.fastq" 
  obigrep --uppercase -l ${params.minLen} -p 'forward_tag is not None and reverse_tag is not None' "${sample_id}_QF_Dmux.fastq" > "${sample_id}_QF_Dmux_minLF.fastq"
  """
}



// Process 4: Split the file based on samples
process split_samples {
  label 'obitools'

  publishDir "04_splitSamples_${sample_id}", mode: params.publish_dir_mode

  input:
    tuple val(sample_id), path(fastqs) 
  
  output:
    tuple val(sample_id), path("$sample_id/*.fastq") 

  script:
  sample_tag = params.illumina_demultiplexed ? "illumina_sample" : "sample"
  """
  mkdir ${sample_id}
  obisplit --uppercase -t ${sample_tag} -u "noSampleID.fastq" $fastqs
  mv *.fastq ${sample_id}
  mv ${sample_id}/$fastqs ..
  mv ${sample_id}/noSampleID.fastq  noSampleID.fastq.ignore
  """
}

/* (vsearch_ch,usearch_ch) = mode == 'vsearch' ? [split_ch, Channel.empty()] : [Channel.empty(), split_ch] */
// Process 5: Relabel file for vsearch
process relabel_vsearch {
  label 'vsearch'

  publishDir "05_${task.process}_${sample_id}", mode: params.publish_dir_mode

  input:
    tuple val(sample_id), path(fastqs) //from vsearch_ch

  output:
    tuple val(sample_id), path("*_upper.fasta"), emit: 'main' //into relabel_ch_vsearch
    /* tuple val(sample_id), path("*.relabeled.fastq"), path("CountOfSeq.txt"), path("*_relabeled4Usearch.fastq"), emit: 'info' //into addition_ch_vsearch */
  
  script:
  """
  for files in ${fastqs}; do
    label=\$(echo \$files | cut -d '/' -f 3 | cut -d '.' -f 1)
    vsearch --threads ${task.cpus} --fastx_filter \$files --relabel \${label}. --fastqout \${label}.relabeled.fastq
  done

  for files in *.relabeled.fastq; do
    name=\$(echo \$files | cut -d '/' -f '2' | cut -d '.' -f 1)
    echo \${name} >> CountOfSeq.txt
    grep "^@\${name}" \$files | wc -l >> CountOfSeq.txt
  done

  cat *.relabeled.fastq > "${sample_id}_QF_Dmux_minLF_relabeled4Usearch.fastq"

  # this function doesn't appear to be supported in vsearch, 
  # but it's not actually used anywher so who cares?
  # ${params.usearch64} -fastx_get_sample_names *_relabeled4Usearch.fastq -output sample.txt

  vsearch --threads ${task.cpus} --fastx_filter *_relabeled4Usearch.fastq --fastaout ${sample_id}.fasta
  # TODO: it's likely that we don't have to do this step if we just pass 'uppercase' to obitools
  awk '/^>/ {print(\$0)}; /^[^>]/ {print(toupper(\$0))}' *.fasta > ${sample_id}_upper.fasta
  """
}

process relabel_usearch {

  label 'usearch'

  publishDir "05_${task.process}_${sample_id}", mode: params.publish_dir_mode

  input:
    tuple val(sample_id), path(fastqs) 

  output:
    tuple val(sample_id), path("*_upper.fasta"), emit: 'main'
    /* tuple val(sample_id), path("*.relabeled.fastq"), path("CountOfSeq.txt"), path("*_relabeled4Usearch.fastq"), emit: 'info'  */

  script:
  if (params.mode == 'usearch32') {
    """
    for files in ${fastqs}; do
      label=\$(echo \$files | cut -d '/' -f 3 | cut -d '.' -f 1)
      usearch -fastx_relabel \$files -prefix \${label}. -fastqout \${label}.relabeled.fastq 
    done

    for files in *.relabeled.fastq; do
      name=\$(echo \$files | cut -d '/' -f '2' | cut -d '.' -f 1)
      echo \${name} >> CountOfSeq.txt
      grep "^@\${name}" \$files | wc -l >> CountOfSeq.txt
    done 

    cat *.relabeled.fastq > "${sample_id}_QF_Dmux_minLF_relabeled4Usearch.fastq"

    usearch -fastx_get_sample_names *_relabeled4Usearch.fastq -output sample.txt

    usearch -fastq_filter *_relabeled4Usearch.fastq -fastaout ${sample_id}.fasta

    awk '/^>/ {print(\$0)}; /^[^>]/ {print(toupper(\$0))}' *.fasta > ${sample_id}_upper.fasta

    """
  }
  else if (params.mode == 'usearch64') {
    """
    for files in ${fastqs}; do
      label=\$(echo \$files | cut -d '/' -f 3 | cut -d '.' -f 1)
      ${params.usearch64} -fastx_relabel \$files -prefix \${label}. -fastqout \${label}.relabeled.fastq
    done

    for files in *.relabeled.fastq; do
      name=\$(echo \$files | cut -d '/' -f '2' | cut -d '.' -f 1)
      echo \${name} >> CountOfSeq.txt
      grep "^@\${name}" \$files | wc -l >> CountOfSeq.txt
    done

    cat *.relabeled.fastq > "${sample_id}_QF_Dmux_minLF_relabeled4Usearch.fastq"

    ${params.usearch64} -fastx_get_sample_names *_relabeled4Usearch.fastq -output sample.txt

    ${params.usearch64} -fastq_filter *_relabeled4Usearch.fastq -fastaout ${sample_id}.fasta

    awk '/^>/ {print(\$0)}; /^[^>]/ {print(toupper(\$0))}' *.fasta > ${sample_id}_upper.fasta
    """
  }
}

/* ch = mode == 'vsearch' ? relabel_ch_vsearch : relabel_ch_usearch */
/* // redirecting channel choice depending on whether or not --skipDemux parameter is set */
/* demux_channel = (params.skipDemux ? name_ch.combine(dmuxed_relabeled_input_ch) : ch) */
  /* (demux_vsearch, demux_usearch) = mode == 'vsearch' ? [demux_channel,Channel.empty()] : [Channel.empty(),demux_channel] */


// Process 6: Dereplication, ZOTUs creation, ZOTU table creation (vsearch version)
process derep_vsearch {
  label 'vsearch'

  publishDir "06_${task.process}_${sample_id}", mode: params.publish_dir_mode

  input:
    tuple val(sample_id), path(upper_fasta) 

  output:
    tuple val(sample_id), path("${sample_id}_Unq.fasta"), path("${sample_id}_zotus.fasta"), path("zotuTable.txt") 

  when:
  !params.onlyDemux

  script:
  """
  vsearch --threads ${task.cpus} --derep_fulllength ${upper_fasta} --sizeout --output "${sample_id}_Unq.fasta"
  vsearch --threads ${task.cpus} --cluster_unoise "${sample_id}_Unq.fasta" --centroids "${sample_id}_centroids.fasta" --minsize ${params.minsize}	   
  vsearch --threads ${task.cpus} --uchime3_denovo "${sample_id}_centroids.fasta" --nonchimeras "${sample_id}_zotus.fasta" --relabel zotu 
  vsearch --threads ${task.cpus} --usearch_global ${upper_fasta} --db "${sample_id}_zotus.fasta" --id 0.97 --otutabout zotuTable.txt
  """
}

process derep_usearch {
  label 'usearch'

  publishDir "06_${task.process}_${sample_id}", mode: params.publish_dir_mode

  input:
    tuple val(sample_id), path(upper_fasta) 

  output:
    tuple val(sample_id), path("${sample_id}_Unq.fasta"), path("${sample_id}_zotus.fasta"), path("zotuTable.txt") 

  when:
    !params.onlyDemux

  script:
  if(params.mode == 'usearch32')
  {
    """
    usearch -fastx_uniques ${upper_fasta} -sizeout -fastaout "${sample_id}_Unq.fasta"
    usearch -unoise3 "${sample_id}_Unq.fasta"  -zotus "${sample_id}_zotus.fasta" -tabbedout "${sample_id}_Unq_unoise3.txt" -minsize ${params.minsize}
    usearch -otutab ${upper_fasta} -zotus ${sample_id}_zotus.fasta -otutabout zotuTable.txt -mapout zmap.txt
    """
  }
  else if(params.mode == 'usearch64')
  {
    """
    ${params.usearch64} -fastx_uniques ${upper_fasta} -sizeout -fastaout "${sample_id}_Unq.fasta"
    ${params.usearch64} -unoise3 "${sample_id}_Unq.fasta"  -zotus "${sample_id}_zotus.fasta" -tabbedout "${sample_id}_Unq_unoise3.txt" -minsize ${params.minsize}
    ${params.usearch64} -otutab ${upper_fasta} -zotus ${sample_id}_zotus.fasta -otutabout zotuTable.txt -mapout zmap.txt
    """
  }
}

// Process 7: Blast
process blast {
  label 'blast'

  cpus { params.max_cpus / 2 }

  publishDir "07_${task.process}_${sample_id}", mode: params.publish_dir_mode

  input:
    tuple val(sample_id), path(a), path(zotus_fasta), path(zotuTable) 

  output:
    tuple val(sample_id), path("${sample_id}_blast_Result.tab"), path(zotuTable), path("match_list.txt") 

  script:
  """
  export BLASTDB="\$(dirname \"${params.blast_db}\")"
  blastn -task ${params.blast_task} \
    -db "${params.blast_db} ${params.custom_db}" \
    -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" \
    -perc_identity ${params.perc_identity} -evalue ${params.evalue} \
    -best_hit_score_edge 0.05 -best_hit_overhang 0.25 \
    -qcov_hsp_perc ${params.qcov} -max_target_seqs ${params.maxTarSeq} \
    -query ${zotus_fasta} -out ${sample_id}_blast_Result.tab \
    -num_threads ${task.cpus}

  # do the blast search we need for LULU
  makeblastdb -in ${zotus_fasta} -parse_seqids -dbtype nucl -out ${sample_id}_zotus
  blastn -db ${sample_id}_zotus \
    -outfmt "6 qseqid sseqid pident" \
    -out match_list.txt -qcov_hsp_perc 80 \
    -perc_identity 84 -query ${zotus_fasta} \
    -num_threads ${task.cpus}

  """
}

// Process 8: LULU curation
process lulu {
  label 'lulu'

  publishDir "08_${task.process}_${sample_id}_minMatch${params.minMatch_lulu}", mode: params.publish_dir_mode

  input:
    tuple val(sample_id), path(a), path(zotuTable), path(match_list) 

  output:
  tuple path("curated_zotuTable.tab"), path("lulu_zotu_map.tab")

  script:
  """
  lulu.R ${params.minMatch_lulu}
  """
}


// Process 9: Taxonomy assignment with LCA
process assign_taxonomy {
  label 'lca_python3'
    publishDir "09_${task.process}_${lca_output}_qCov${lca_qcov}_id${lca_pid}_diff${lca_diff}", mode: params.publish_dir_mode

    input:
    tuple path(table), path(blastRes) from zotuTable_ch.combine(blastFile_ch)
    file lcaScript from lca_script

    output:
    tuple path("interMediate_res.tab"), path("${lca_output}_qCov${lca_qcov}_id${lca_pid}_diff${lca_diff}.tab") into last_ch


    script:
    """
    python3 $lca_script ${table} ${blastRes} ${lca_qcov} ${lca_pid} ${lca_diff} ${lca_output}_qCov${lca_qcov}_id${lca_pid}_diff${lca_diff}.tab
    """
}  

// for some reason we have to do this in here also, or else
// the case conversion business doesn't seem to work
params.r1 = "R1"
params.r2 = "R2"
params.fwd = ""
params.rev = ""
params.baseDir = "."
params.reads   = "*_{R1,R2}.fastq"
params.barcode = "*_bc.txt"
params.minQuality = "20"
params.minAlignLeng = "12"
params.minLen = "50"
params.primer_mismatch = "2"
params.minsize = "8"
params.maxTarSeq = "10"
params.perc_identity = "95"
params.evalue = "1e-3"
params.qcov = "100"
params.lulu = "lulu.R"
params.mode = "usearch32"  
params.usearch64 = ""   
params.blast_db = ""
params.custom_db = ""
params.blast_task = "blastn"
params.publish_dir_mode = "symlink" 
params.bindDir = ""
params.singularityDir = ""
params.prefix = "seq"
params.illumina_demultiplexed = false
params.concat = false
params.remove_ambiguous_tags = false
params.onlyDemux = false
params.skipDemux = false
params.demuxedInput = "*.fasta"
params.skipFastqc = false
params.taxonomyAssignment = false
params.lca_script = "LCA_taxonomyAssignment_scripts/runAssign_collapsedTaxonomy.py"
params.zotuTable = ""
params.blastFile = ""
params.lca_qcov = "100"
params.lca_pid = "97"
params.lca_diff = "1"
params.lca_output = "lca_taxAssigned_results"
params.minMatch_lulu="84"
params.test = false
params.help = false
params.max_memory = Runtime.runtime.maxMemory()
params.max_cpus = Runtime.runtime.availableProcessors()
params.max_time = "240.h"

workflow {
  // show help message and bail
  if (params.help) {
    usage()
    exit(0)
  }

  if (!params.skipDemux) {
    if (params.fwd != "" && params.rev != "") {
      pattern = "${params.baseDir}/{${params.fwd},${params.rev}}/*{${params.r1},${params.r2}}*.fastq"
    } else {
      pattern = params.reads
    }

    // load the reads and make sure they're in the right order because the
    // files will be sorted alphabetically this is an extremely niche issue
    // because file are basically always R1/R2 but what if they're
    // forwards/backwards or something like that?

    // I'm not even certain the order matters, but I think it does because we 
    // send them to --file1 and --file2 of AdapterRemoval
    reads = Channel.fromFilePairs(pattern).map { key,f -> [key,f[0] =~ /${params.r1}/ ? [f[0],f[1]] : [f[1],f[0]]]  }
    barcodes = Channel.fromPath(params.barcode, checkIfExists: true)

    // if the sequences are already demultiplexed by illumina, we'll
    // annotate them with their sample names and optionally filter out
    // sequences with ambiguous tags. then we recombine them into one forward
    // and one reverse file and continue with those as our reads
    if (params.illumina_demultiplexed) {
      // this is a whole bunch of files potentially
      reads_annotated = annotate_demultiplexed(reads)

      // here's where we smash all the annotated reads into one file for each direction
      fwd = reads_annotated.flatMap {k, f1, f2 -> f1}.collectFile(name: 'fwd.fastq',newLine: false, sort: 'index')
      rev = reads_annotated.flatMap {k, f1, f2 -> f2}.collectFile(name: 'rev.fastq',newLine: false, sort: 'index')
      reads = Channel.of(params.prefix).
        combine(fwd).
        combine(rev).
        map { k,f1,f2 -> [k,[f1,f2]] }
      // at this point we should only have two files
    }

    // do initial fastqc step
    if (!params.skipFastqc) {
      initial_fastqc(reads,'initial')
    }

    // do quality filtering and/or paired-end merge
    reads_filtered_merged = filter_merge(reads)

    // do after-filtering fastqc step
    if (!params.skipFastqc) {
      merged_fastqc(reads_filtered_merged,'filtered_trimmed')
    }

    // demultiplex samples (if necessary)
    demultiplexed = demultiplex(reads_filtered_merged.combine(barcodes)).groupTuple()

    // do length filtering
    length_filtered = filter_length(demultiplexed)

    // split reads by sample
    samples_split = split_samples(length_filtered)

    /* to_dereplicate = params.mode == 'vsearch' ?  */
      /* relabel_vsearch(samples_split) : */
      /* relabel_usearch(samples_split) */
    if (params.mode == 'vsearch') {
      to_dereplicate = relabel_vsearch(samples_split)
    } else {
      to_dereplicate = relabel_usearch(samples_split)
    }
  } else {
    to_dereplicate = Channel.value(params.demuxedInput).
      combine(Channel.fromPath(params.demuxedInput))
  }

  dereplicated = params.mode == 'vsearch' ? 
    derep_vsearch(to_dereplicate) :
    derep_usearch(to_dereplicate)
  /* if (params.mode == 'vsearch') { */
    /* dereplicated = derep_vsearch(to_dereplicate) */
  /* } else { */
    /* dereplicated = derep_usearch(to_dereplicate) */
  /* } */

  blast_results = blast(dereplicated)
  lulu_results = lulu(blast_results)
}
