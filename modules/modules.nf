// run fastqc process
process fastqc {
  label 'fastqc'
  label 'process_single'

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
  label 'process_single'

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