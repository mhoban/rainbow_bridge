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

// get a file from a URL
process get_web {
  label 'process_single'

  input:
    tuple val(location), path(localfile)
  output:
    path(localfile)

  exec:
    def classifier = task.workDir / localfile 
    def url = new URL(location)
    url.withInputStream { stream -> classifier << stream }
}

// extract arbitrary files from a zip archive
process extract_zip {
  container 'quay.io/biocontainers/unzip:6.0'
  label 'process_single'

  input:
    tuple path(zipfile), val(f)
  output:
    path(f)
  
  script:
  """
  unzip -p ${zipfile} ${f} > ${f}
  """
}

// extract files from a .tar.gz archive
process extract_targz {
  label 'shell'
  label 'process_single'

  input:
    tuple path(archive), val(to_extract)
  
  output:
    path(to_extract)

  script:
  """
  gunzip -c ${archive} | tar x ${to_extract}
  """
}