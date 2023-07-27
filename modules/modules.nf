// assign/collapse taxonomy using the R port of the original python script
process r_taxonomy {
  label 'phyloseq'

  publishDir { 
    p = params.insect ? "a" : ""
    return params.illuminaDemultiplexed ? "09${p}_collapsed_taxonomy_r" : "10${p}_collapsed_taxonomy_r" 
  }, mode: params.publishMode

  input:
    tuple path(zotu_table), path(blast_result), path(lineage), val(name)

  output:
    tuple path("${name}_intermediate_table.tab"), path("${name}_taxonomy_collapsed.tab")


  script:
  pf = params.filterUncultured ? "-f" : ""
  """
  collapse_taxonomy.R \
    --qcov ${params.lcaQcov} \
    --pid ${params.lcaPid} \
    --diff ${params.lcaDiff} \
    ${pf} \
    --intermediate "${name}_intermediate_table.tab" \
    ${blast_result} ${zotu_table} ${lineage} "${name}_taxonomy_collapsed.tab"
  """
}  

// assign/collapse taxonomy using the original python script
process taxonomy {
  label 'python3'

  publishDir { 
    p = params.insect ? "a" : ""
    return params.illuminaDemultiplexed ? "09${p}_collapsed_taxonomy" : "10${p}_collapsed_taxonomy" 
  }, mode: params.publishMode

  input:
    tuple path(zotu_table), path(blast_result), val(name)

  output:
    tuple path("${name}_intermediate_table.tab"), path("${name}_taxonomy_collapsed.tab")


  script:
  """
  runAssign_collapsedTaxonomy.py ${zotu_table} ${blast_result} ${params.lcaQcov} ${params.lcaPid} ${params.lcaDiff} ${name}_taxonomy_collapsed.tab
  mv interMediate_res.tab ${name}_intermediate_table.tab
  """
}  

// run fastqc process
process fastqc {
  label 'fastqc'

  publishDir { "00_fastqc_${step}" }, mode: params.publishMode

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

  publishDir { "00_fastqc_${step}" }, mode: params.publishMode

  input:
    tuple path(fastqc_files), val(step)
  output: 
    tuple path('multiqc_report.html'), path('multiqc_data/*')

  script:
  """
  multiqc .
  """
}
