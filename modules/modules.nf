// assign/collapse taxonomy using the R port of the original python script
process lca {
  label 'r'

  publishDir "${params.outDir}/taxonomy/lca/qcov${params.lcaQcov}_pid${params.lcaPid}_diff${params.lcaDiff}", mode: params.publishMode 

  input:
    tuple path(blast_result), path(lineage), path(merged), val(curated)

  output:
    tuple path("lca_intermediate*.tsv"), path("lca_taxonomy*.tsv"), emit: result
    path 'lca_settings.txt'


  script:
  pf = params.keepUncultured ? "-u" : ""
  c = curated ? "_lulu" : ""
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
    --merged ${merged} \
    --dropped ${params.dropped} \
    ${pf} \
    --intermediate "lca_intermediate${c}.tsv" \
    ${blast_result} ${lineage} "lca_taxonomy${c}.tsv"
  """
}  

// run fastqc process
process fastqc {
  label 'fastqc'

  publishDir "${params.outDir}/fastqc/${step}", mode: params.publishMode

  input:
    tuple val(step), val(sample_id), path(read) 

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
