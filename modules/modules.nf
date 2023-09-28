// assign/collapse taxonomy using the R port of the original python script
process r_taxonomy {
  label 'phyloseq'

  publishDir { 
    p = params.insect ? "a" : ""
    if (!standalone) {
      // number output directory if it's part of the pipeline
      dir =  params.illuminaDemultiplexed ? "09${p}_taxonomy" : "10${p}_taxonomy" 
    } else {
      // otherwise we're standalone, so don't
      dir = "taxonomy"
    }
    // include taxonomy parameters in output directory name
    return "${dir}_q${params.lcaQcov}_p${params.lcaPid}_d${params.lcaDiff}"
  }, mode: params.publishMode

  input:
    tuple path(zotu_table), path(blast_result), path(lineage), val(name), val(standalone)

  output:
    tuple path("${name}_intermediate_r.tab"), path("${name}_taxonomy_r.tab")


  script:
  pf = params.filterUncultured ? "-f" : ""
  """
  collapse_taxonomy.R \
    --qcov ${params.lcaQcov} \
    --pid ${params.lcaPid} \
    --diff ${params.lcaDiff} \
    ${pf} \
    --intermediate "${name}_intermediate_r.tab" \
    ${blast_result} ${zotu_table} ${lineage} "${name}_taxonomy_r.tab"

  """
}  

// assign/collapse taxonomy using the original python script
process py_taxonomy {
  label 'python3'

  publishDir { 
    p = params.insect ? "a" : ""
    if (!standalone) {
      // number output directory if it's part of the pipeline
      dir =  params.illuminaDemultiplexed ? "09${p}_taxonomy" : "10${p}_taxonomy" 
    } else {
      // otherwise we're standalone, so don't
      dir = "taxonomy"
    }
    // include taxonomy parameters in output directory name
    return "${dir}_q${params.lcaQcov}_p${params.lcaPid}_d${params.lcaDiff}"
  }, mode: params.publishMode

  input:
    tuple path(zotu_table), path(blast_result), path(lineage), val(name), val(standalone)

  output:
    tuple path("${name}_intermediate_py.tab"), path("${name}_taxonomy_py.tab")


  script:
  """
  runAssign_collapsedTaxonomy.py ${zotu_table} ${blast_result} ${params.lcaQcov} ${params.lcaPid} ${params.lcaDiff} ${lineage} ${name}_taxonomy_py.tab 
  mv interMediate_res.tab ${name}_intermediate_py.tab
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
