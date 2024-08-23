import nextflow.Nextflow

class helper {
  static public Map insect_classifiers = [
    'mifish': 'https://www.dropbox.com/s/fv3dpvws6zjvtib/classifier.rds?dl=1',
    'crust16s': 'https://www.dropbox.com/s/9vl9gj3frw7ng1m/classifier.rds?dl=1',
    'fish16s': 'https://www.dropbox.com/s/fvfrd46exdah037/classifier.rds?dl=1',
    '18suni': 'https://www.dropbox.com/s/rmhh1g73jtipagu/classifier.rds?dl=1',
    '18sv4': 'https://www.dropbox.com/s/s315gxuo4p24kx8/classifier.rds?dl=1',
    'p23s': 'https://www.dropbox.com/s/6o8cauqrlgnmwp5/classifier.rds?dl=1',
    'mlcoiint': 'https://www.dropbox.com/s/dvnrhnfmo727774/classifier.rds?dl=1',
    'scl5.8S': 'https://www.dropbox.com/s/f07cka6308ebk2o/classifier.rds?dl=1'
  ]

  // get common characters from left side of two strings
  static public String common(one, two) {
    def comm = ""
    for (def i=0; i< one.size(); i++) {
      if (one[i] == two[i]) comm += one[i]
      else break
    }
    return comm
  }

	static public String basename(f) {
		return (new File(f)).getBaseName()
	}

	static public boolean file_exists(f) {
		return (new File(f)).exists()
	}

	static public boolean is_dir(d) {
		return (new File(d)).isDirectory()
	}

	static public boolean executable(f) {
		return (new File(f)).canExecute()
	}

	static public String get_env(n) {
		
		return System.getenv(n) ? System.getenv(n) : ""
	}
  
  static public boolean is_list(object) {    
      [Collection, Object[]].any { it.isAssignableFrom(object.getClass())  }
  }

	static public void demuxed_example() {
		System.out.println("""
		>sample1.1
		AGCGTCCGATGACTGACTGACTAGCT
		>sample1.2
		TACGTACGATCGACGAGTCTACGACTACTGAC
		>sample1.3
		TGACTGATCGTACTATCAGAGCTATCATCGACTATCATCGATC
		>sample2.1
		ATCGTACTACTAGCGACGAGTCATCACGACGTACTAGTCGA
		>sample2.2
		CATGCGACGTACGTACTATCATCATCGAGCAGCTATATATCGATGGTACTAGCTGAC
		>sample2.3
		TGACTGATCGTACTATCAGAGCTATCATCGACTATCATCGATC
		>sample3.1
		AGCGTCCGATGACTGACTGACTAGCT
		>sample3.2
		ATCGTACTACTAGCGACGAGTCATCACGACGTACTAGTCGA
		>sample3.3
		CATGCGACGTACGTACTATCATCATCGAGCAGCTATATATCGATGGTACTAGCTGAC
		""".stripIndent())
	}

  static public void usage(params) {
		System.out.println("""
    Usage: rainbow_bridge.nf [options]

    General options:
      --project [project]          Project name, applied to various output filenames (default: ${params.project}) 
      --barcode [file]             (required) Barcode file. Format must match OBITools requirements
                                   (see https://pythonhosted.org/OBITools/scripts/ngsfilter.html)
                                   To denote multiple barcode files, this may be a glob, but it must
                                   be surrounded by quotations (e.g. 'barcode*.tab'). 
                                   For previously-demultiplexed datasets, use ':' for barcode pairs
                                   Primer sequences are still used for primer-mismatch comparisons
      --publish-mode [mode]        Specify how nextflow places files in output directories 
                                   (default: symlink)
      --fastqc                     Output FastQC reports for pre and post filter/merge steps 
      --split                      Split input fastq files and process in parallel
                                   (not compatible with --illumina-demultiplexed)
      --split-by                   Number of sequences per split fastq chunk (default: ${params.splitBy})
      --sample-map [mapfile]       (Optional) A headerless tab-separated file mapping sample names to sequence-read
                                   filenames. Paired-end runs include both forward and reverse reads. Example map:
                                   ---
                                   sample1	B1_S7_L001_R1_001.fastq	B1_S7_L001_R2_001.fastq
                                   sample2	B2_S8_L001_R1_001.fastq	B2_S8_L001_R2_001.fastq
                                   sample3	CL1_S2_L001_R1_001.fastq	CL1_S2_L001_R2_001.fastq
                                   sample4	CL2_S3_L001_R1_001.fastq	CL2_S3_L001_R2_001.fastq 
                                   ---
                                   NOTE: if your fastq files are gzipped, DO NOT include the .gz extension 
                                   in your sample map file, because the files will be unzipped 
                                   (and .gz extension stripped) BEFORE sample IDs are remapped

    For single-end sequencing runs:
      --single                     Specify single-ended sequencing run (required)
      --reads [file/glob]          Location of sequencing read(s). For runs that have already been 
                                   demultiplexed by the sequencer, this may be a glob (e.g., '*.fastq')

    For paired-end sequencing runs:  
      --paired                     Specify paired-end sequencing run (required)

    To specify the location of paired-end sequence files, the following methods are availble

    Resovle paird-end reads locations directly using globs:
      --reads [glob]               If --reads is a glob, attempt to resolve paired-end reads directly
                                   e.g., '/data/run1/*{R1,R2}*.fastq.gz'
      --fwd [glob], --rev [glob]   Resolve forward and reverse reads directly using globs
                                   e.g., --fwd 'r1/*R1*.fastq' --rev 'r2/*R2*.fastq'

    Specify location of paired-end reads with directories, using the following patterns:
      <reads>/*<r1>|<r2>*.f*q*
      <reads>/<fwd>|<rev>/*<r1>|<r2>*.f*q*
      <fwd>|<rev>/*<r1>|<r2>*.f*q*

      --reads [dir]                 Location (directory) where sequence reads can be found (default: .) 
      --fwd [dir]                   (Optional) forward reads directory (default: ${params.fwd})
                                    For runs that have NOT been demultiplexed, --fwd may also point
                                    directly at the raw forward reads (R1) sequence file.
      --rev [dir]                   (Optional) reverse reads directory (default: ${params.rev})
                                    For runs that have NOT been demultiplexed, --rev may also point
                                    directly at the raw reverse reads (R2) sequence file.
      --r1 [pattern]                Pattern distinguishing forward read files (default: ${params.r1})
      --r2 [pattern]                Pattern distinguishing reverse read files (default: ${params.r2})

    Length, quality, and merge settings:
      --min-quality [num]           Minimum Phred score for sequence retention 
                                    (default: ${params.minQuality})
      --min-align-len [num]         Minimum sequence overlap when merging forward/reverse reads
                                    (default: ${params.minAlignLen})
      --min-len [num]               Minimum overall sequence length (default: ${params.minLen})

    BLAST settings (one or more database is required, unless skipping):
      --skip-blast                  Skip BLAST searches
      --blast-db [blastdb]          Location of BLAST database (path AND name).
                                    e.g., '/drive1/blast/custom_db', where the database files are named
                                    things like custom_db.ndb, custom_db.nhr, custom_db.nin, etc.
                                    If the \$FLOW_BLAST environment variable points to a BLAST database, 
                                    the pipeline will include it in the list of databases to search.
                                    To specify multiple BLAST databases, pass them as a list in the
                                    parameters file (see README for more information).
      --blast-taxdb [dir]           Specify the location of the taxdb.* files, for BLAST taxon name assignment.
                                    If unspecified or missing from the --blast-db directory, the pipeline 
                                    will download these files from the NCBI server.
      --ignore-blast-env            Ignore the value of the \$FLOW_BLAST environment variable
      --blast-task [task]           Set blast+ task (default: "blastn")
      --max-query-results [num]     Maxmimum number of BLAST results to return per zOTU (default: 10)
      --percent-identity [num]      Minimum percent identity of matches to report (default: 95)
      --evalue [num]                Expectation value threshold for saving hits (default: 0.001)
      --qcov [num]                  Percent query coverage per hsp (default: 100)

    LCA taxonomy collapse:
      --collapse-taxonomy           Collapse assigned BLAST results by least common ancestor (LCA)
      --standalone-taxonomy         Run LCA script as standalone
      --blast-file [file]           Blast result table (only for standalone LCA assignment)
      --zotu-table [file]           zOTU table file (only for standalone LCA assignment)
      --taxdump [file]              Previously downloaded NCBI taxonomy dump zip archive (new_taxdump) (leave blank to download)
      --lca-qcov [num]              Minimum query coverage for LCA taxonomy assignment (default: 100)
      --lca-pid [num]               Minimum percent identity for LCA taxonomy assignment (default: 97)
      --lca-diff [num]              The difference between percent identities (when query coverage is
                                    identical) where species-level taxonomy is retained (default: 1)
      --keep-uncultured             Keep sequences that are listed as 'uncultured', 'environmental sample',
                                    'synthetic', or 'clone'
      --dropped [str]               Placeholder text for dropped taxonomic levels (use "NA" for blank/NA)

    Insect taxonomy classification:
      --insect [classifier]         Perform taxonomy assignment using insect
                                    Accepted values of [classifier] are:
                                    - Filename to local .rds containing classifier model
                                    - One of the following (case-insensitive) primer names: 
                                      MiFish, Crust16S, Fish16S, 18SUni, 18SV4, p23S, mlCOIint, SCL5.8S
                                      (see https://github.com/shaunpwilkinson/insect#classifying-sequences)
      --insect-threshold [num]      Minimum Akaike weight for the recursive classification procedure to continue 
                                    toward the leaves of the tree (default: ${params.insectThreshold})
      --insect-offset [num]         Log-odds score offset parameter governing whether the minimum score is 
                                    met at each node (default: ${params.insectOffset})
      --insect-min-count [num]      Minimum number of training sequences belonging to a selected child node 
                                    for the classification to progress (default: ${params.insectMinCount})
      --insect-ping [num]           Numeric (between 0 and 1) indicating whether a nearest neighbor search should 
                                    be carried out, and if so, what the minimum distance to the nearest neighbor 
                                    should be for the the recursive classification algorithm to be skipped (default: ${params.insectPing})

    Generating phyloseq objects:
      --phyloseq                    Create phyloseq object (requires --collapse-taxonomy)
      --metadata [file]             Comma or tab-separated sample metadata file (required)
      --taxonomy [tax]                Taxonomic classifaction scheme. May be pipeline-generated or user supplied
                                    (acceptable options: lca, insect, <filename>; default: ${params.taxonomy})
      --no-tree                     Do not include a phylogenetic tree
      --optimize-tree               Attempt to optimize generation of the tree (may take a long time)

    Data cleanup & rarefaction:
      --abundance-filter            Filter reads by minimum relative abundance
      --abundance-threshold [num]   Minimum relative abundance below which counts are set to zero
      --filter-minimum              Filter samples by minimum read count
      --min-reads [num]             Total read count below which samples will be filtered
      --rarefy                      Rarefy samples to minimum depth
      --rarefaction-method [method] Rarefaction method (available options: perm, phyloseq)
      --permutations [num]          For 'perm' rarefaction method, number of permutations

    Decontamination & taxon filtering
      --taxon-remap [file]          Taxonomy remap file (.csv or .tsv)
      --taxon-filter [file]         Taxonomy filter file (.csv or .tsv)
      --taxon-priority [str]        Priority when taxonomies disagree (values: lca or insect)
      --controls [file]             File listing negative control sample IDs
      --control-action [action]     Negative control decontamination method ('remove', 'subtract', or 'decontam')
      --control-threshold [num]     Minimum read threshold or decontam threshold value
      --decontam-method [method]    Method passed to isContaminant function of decontam (default: 'auto')
      --dna-concentration [file]    File specifying DNA concentrations by sample ID

    Demultiplexing and sequence matching:
      --illumina-demultiplexed      Sequencing run has already been demultiplexed by the sequencer

      --remove-ambiguous-indices    Removes reads with ambiguous indices in the header (i.e., not A,G,C,T)
                                    (only applies to previously-demultiplexed runs with indices in header)
      --demuxed-fasta [file]        Skip demultiplexing step and use supplied FASTA 
                                    (must be in usearch/vsearch format)
      --demuxed-example             Spit out example usearch/vsearch demultiplexed FASTA format
      --demux-only                  Stop after demultiplexing and splitting raw reads
      --primer-mismatch             Allowed number of mismatched primer bases 
                                    (default: ${params.primerMismatch})
      --skip-primer-match           Skip primer matching (ngsfilter) altogether. 
                                    Use with demultiplexed runs lacking primer sequences.

    Denoising and zOTU inference:  
      --denoiser [tool/path]        Sets the tool used for denoising & chimera removal
                                    accepted options: usearch/usearch32 (equivalent), vsearch, 
                                    path to usearch64 executable (default: vsearch)
      --alpha [num]                 Sets the alpha parameter for the UNOISE3 algorithm (default: ${params.alpha})
      --min-abundance [num]         Minimum sequence abundance for zOTU determination; sequences below threshold will be discarded
                                    (default: ${params.minAbundance}) 
      --usearch                     shortcut for --denoiser usearch

    LULU zOTU curation:
      --skip-lulu                   Skip LULU curation
      --lulu-min-ratio-type [num]   LULU minimum ratio type (accepted values: 'min', 'avg', default: ${params.luluMinRatioType})
      --lulu-min-ratio [num]        LULU minimum ratio (default: ${params.luluMinRatio})
      --lulu-min-match [num]        LULU minimum threshold of sequence similarity to consider zOTUs as spurious (default: ${params.luluMinMatch})
                                    Choose higher values when using markers with lower genetic variation 
                                    and/or few expected PCR and sequencing errors. (default: ${params.luluMinMatch})
      --lulu-min-rc [num]           LULU minimum relative co-occurence rate (default: ${params.luluMinRc})

    Resource allocation:
      --max-memory [mem]            Maximum memory available to nextflow processes, e.g., '8.GB' (default: ${params.maxMemory})
      --max-cpus [num]              Maximum cores available to nextflow processes default: ${params.maxCpus})
      --max-time [time]             Maximum time allocated to each pipeline process, e.g., '2.h' (default: ${params.maxTime})

    Singularity options:
      --bind-dir [dir]              Space-separated list of directories to bind within singularity images
                                    (must be surrounded by quotations if more than one directory)
                                    Note: singularity will attempt to auto-bind all provided host paths
                                    so this option may not be necessary, but try it if you're getting 
                                    weird "file not found" types of errors
      --singularity-cache [dir]     Location to store singularity images. May also be specified
                                    with the environment variable \$NXF_SINGULARITY_CACHEDIR.
                                    (current value: ${get_env("NXF_SINGULARITY_CACHEDIR")})
		""".stripIndent())
  }
}
