// assign/collapse taxonomy using the R port of the original python script
process lca {
  label 'r'

  publishDir "${params.outDir}/taxonomy/lca/qcov${params.lcaQcov}_pid${params.lcaPid}_diff${params.lcaDiff}", mode: params.publishMode 

  input:
    tuple path(blast_result), path('*')

  output:
    tuple path("lca_intermediate*.tsv"), path("lca_taxonomy*.tsv"), emit: result
    path 'lca_settings.txt'


  script:
  pf = params.keepUncultured ? "-u" : ""
  """
  # save settings
  echo "Minimum query coverage %: ${params.lcaQcov}" > lca_settings.txt
  echo "Minimum percent identity: ${params.lcaPid}" >> lca_settings.txt
  echo "Minium percent identity difference: ${params.lcaDiff}" >> lca_settings.txt
  echo "Retain \"uncultured\" results: ${params.keepUncultured ? 'yes' : 'no'}" >> lca_settings.txt

  collapse_taxonomy.R \
    --qcov ${params.lcaQcov} \
    --pid ${params.lcaPid} \
    --diff ${params.lcaDiff} \
    --merged merged.dmp \
    --dropped ${params.dropped} \
    ${pf} \
    --intermediate "lca_intermediate.tsv" \
    ${blast_result} rankedlineage.dmp nodes.dmp "lca_taxonomy.tsv"
  """
}  

// run fastqc process
process fastqc {
  label 'fastqc'

  publishDir "${params.outDir}/fastqc/${step}", mode: params.publishMode

  input:
    tuple val(step), val(key), path(read) 

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
process multiqc {
  label 'multiqc'

  publishDir "${params.outDir}/fastqc/${step}", mode: params.publishMode

  input:
    tuple path(fastqc_files), val(step)
  output: 
    tuple path('multiqc_report.html'), path('multiqc_data/*')

  script:
  """
  multiqc .
  """
}

// get a file from a URL
process get_web {
  input:
    tuple val(location), path(localfile)
  output:
    path(localfile)

  exec:
    classifier = task.workDir / localfile 
    url = new URL(location)
    url.withInputStream { stream -> classifier << stream }
}