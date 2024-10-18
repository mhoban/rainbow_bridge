// assign/collapse taxonomy using the R port of the original python script
process lca {
  label 'r'

  publishDir {
    td = params.standaloneTaxonomy ? 'standalone_lca' : 'lca'
    "${params.outDir}/taxonomy/${td}/qcov${params.lcaQcov}_pid${params.lcaPid}_diff${params.lcaDiff}"
  }, mode: params.publishMode

  input:
    tuple path(blast_result), path(lineage), path('*')

  output:
    tuple path("lca_intermediate*.tsv"), path("lca_taxonomy*.tsv"), emit: result
    path 'lca_settings.txt'


  script:
  pf = []
  params.lcaFilterMaxQcov && pf << "--filter-max-qcov"
  params.lcaCaseInsensitive && pf << "--case-insensitive" 
  """
  # save settings
  echo "Minimum query coverage %: ${params.lcaQcov}" > lca_settings.txt
  echo "Minimum percent identity: ${params.lcaPid}" >> lca_settings.txt
  echo "Minium percent identity difference: ${params.lcaDiff}" >> lca_settings.txt
  echo "Filter to maximum query coverage: ${params.lcaFilterMaxQcov ? 'yes' : 'no'}" >> lca_settings.txt
  echo "Filter taxa by regex: ${params.lcaTaxonFilter}" >> lca_settings.txt
  echo "Taxon filter case sensitive: ${!params.lcaCaseInsensitive ? 'yes' : 'no'}" >> lca_settings.txt

  collapse_taxonomy.R \
    --qcov ${params.lcaQcov} \
    --pid ${params.lcaPid} \
    --diff ${params.lcaDiff} \
    --merged merged.dmp \
    --nodes nodes.dmp \
    --taxid-lineage taxidlineage.dmp \
    --output lca_taxonomy.tsv \
    --dropped "${params.dropped}" \
    --intermediate "lca_intermediate.tsv" \
    ${pf.join(" ")} \
    --taxon-filter "${params.lcaTaxonFilter}" \
    ${blast_result} ${lineage}
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

// extract arbitrary files from a zip archive
process extract_zip {
  label 'shell'

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

  input:
    tuple path(archive), val(to_extract)
  
  output:
    path(to_extract)

  script:
  """
  gunzip -c ${archive} | tar x ${to_extract}
  """
}