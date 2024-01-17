// assign/collapse taxonomy using the R port of the original python script
process r_taxonomy {
  label 'phyloseq'

  publishDir "${params.outDir}/taxonomy/lca", mode: params.publishMode 

  input:
    tuple path(zotu_table), path(blast_result), path(lineage), path(merged), val(curated)

  output:
    tuple path("lca_intermediate*.tab"), path("lca_taxonomy*.tab")


  script:
  pf = params.filterUncultured ? "-f" : ""
  c = curated ? "lulu_" : ""
  """
  collapse_taxonomy.R \
    --qcov ${params.lcaQcov} \
    --pid ${params.lcaPid} \
    --diff ${params.lcaDiff} \
    --merged ${merged} \
    --dropped ${params.dropped} \
    ${pf} \
    --intermediate "lca_intermediate_q${params.lcaQcov}_p${params.lcaPid}_d${params.lcaDiff}_${c}r.tab" \
    ${blast_result} ${zotu_table} ${lineage} "lca_taxonomy_q${params.lcaQcov}_p${params.lcaPid}_d${params.lcaDiff}_${c}r.tab"

  """
}  

// assign/collapse taxonomy using the original python script
process py_taxonomy {
  label 'python3'

  publishDir "${params.outDir}/taxonomy/lca", mode: params.publishMode 

  input:
    tuple path(zotu_table), path(blast_result), path(lineage), path(merged), val(curated)

  output:
    tuple path("lca_intermediate*.tab"), path("lca_taxonomy*.tab")


  script:
  c = curated ? "lulu_" : ""
  """
  runAssign_collapsedTaxonomy.py ${zotu_table} ${blast_result} \
    ${params.lcaQcov} ${params.lcaPid} ${params.lcaDiff} ${lineage} \
    "lca_taxonomy_q${params.lcaQcov}_p${params.lcaPid}_d${params.lcaDiff}_${c}py.tab" 
  mv interMediate_res.tab "lca_intermediate_q${params.lcaQcov}_p${params.lcaPid}_d${params.lcaDiff}_${c}py.tab"
  """
}  

// run fastqc process
process fastqc {
  label 'fastqc'

  publishDir "${params.outDir}/fastqc/fastqc_${step}", mode: params.publishMode

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

  publishDir "${params.outDir}/fastqc/fastqc_${step}", mode: params.publishMode

  input:
    tuple path(fastqc_files), val(step)
  output: 
    tuple path('multiqc_report.html'), path('multiqc_data/*')

  script:
  """
  multiqc .
  """
}
