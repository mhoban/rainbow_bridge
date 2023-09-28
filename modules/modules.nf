// assign/collapse taxonomy using the R port of the original python script
process r_taxonomy {
  label 'phyloseq'

  publishDir { 
    p = params.insect ? "a" : ""
    if (!params.standaloneTaxonomy) {
      num = (params.illuminaDemultiplexed ? 8 : 9) + (params.skipLulu ? 0 : 1)
      num = String.format("%02d",num)

      // number output directory if it's part of the pipeline
      dir =  params.illuminaDemultiplexed ? "${num}${p}_taxonomy" : "${num}${p}_taxonomy" 
    } else {
      // otherwise we're standalone, so don't
      dir = "taxonomy"
    }
    return dir
  }, mode: params.publishMode

  input:
    tuple path(zotu_table), path(blast_result), path(lineage), path(merged), val(curated)

  output:
    tuple path("intermediate*.tab"), path("taxonomy*.tab")


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
    --intermediate "intermediate_q${params.lcaQcov}_p${params.lcaPid}_d${params.lcaDiff}_${c}r.tab" \
    ${blast_result} ${zotu_table} ${lineage} "taxonomy_q${params.lcaQcov}_p${params.lcaPid}_d${params.lcaDiff}_${c}r.tab"

  """
}  

// assign/collapse taxonomy using the original python script
process py_taxonomy {
  label 'python3'

  publishDir { 
    p = params.insect ? "a" : ""
    if (!params.standaloneTaxonomy) {
      // number output directory if it's part of the pipeline
      num = (params.illuminaDemultiplexed ? 8 : 9) + (params.skipLulu ? 0 : 1)
      num = String.format("%02d",num)

      // number output directory if it's part of the pipeline
      dir =  params.illuminaDemultiplexed ? "${num}${p}_taxonomy" : "${num}${p}_taxonomy" 
    } else {
      // otherwise we're standalone, so don't
      dir = "taxonomy"
    }
    return dir
  }, mode: params.publishMode

  input:
    tuple path(zotu_table), path(blast_result), path(lineage), path(merged), val(curated)

  output:
    tuple path("intermediate*.tab"), path("taxonomy*.tab")


  script:
  c = curated ? "lulu_" : ""
  """
  runAssign_collapsedTaxonomy.py ${zotu_table} ${blast_result} ${params.lcaQcov} ${params.lcaPid} ${params.lcaDiff} ${lineage} "taxonomy_q${params.lcaQcov}_p${params.lcaPid}_d${params.lcaDiff}_${c}py.tab" 
  mv interMediate_res.tab "intermediate_q${params.lcaQcov}_p${params.lcaPid}_d${params.lcaDiff}_${c}py.tab"
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
