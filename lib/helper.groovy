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
    'scl5.8s': 'https://www.dropbox.com/s/f07cka6308ebk2o/classifier.rds?dl=1'
  ]

  // get common characters from left side of two strings
  static public String common(one, two) {
    def comm = ""
    for (def i=0; i < Math.min(one.size(),two.size()); i++) {
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

  static public boolean is_url(url) throws MalformedURLException, URISyntaxException {
    try {
      new URL(url).toURI();
      return true;
    } catch (MalformedURLException e) {
      return false;
    } catch (URISyntaxException e) {
      return false;
    }
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

  static public void usage(params) {
		System.out.println("""
    Usage: rainbow_bridge.nf [options]

    Useful nextflow options (note single leading dash):
      -params-file [file]          Load rainbow_bridge options from YAML or json parameters file
      -N [email address]           Notify by email on pipeline completion or error

    General options:
      --demultiplexed-by [option]  (required) Specify demultiplexing strategy. Accepted options are 
                                   `index`, `barcode`, `combined` (default: ${params.demultiplexedBy})
      --barcode [file]             (required) Barcode file. Format must match OBITools requirements
                                   (see https://pythonhosted.org/OBITools/scripts/ngsfilter.html)
                                   To denote multiple barcode files, this may be a glob, but it must
                                   be surrounded by quotations (e.g. 'barcode*.tab'). 
                                   For previously-demultiplexed datasets, use ':' for barcode pairs
                                   Primer sequences are still used for primer-mismatch comparisons
                                   For `combined` demultiplexing strategy, first column of barcode file
                                   must match the underscore-delimited prefix of your sequence read files.
                                   (See README for more details).
      --fwd-primer [seq]           Forward PCR primer (to trim)
      --reverse-primer [seq]       Reverse PCR primer (to trim)
      --free-primers               Primer sequences are unanchored when using cutadapt 
      --project [project]          Project name, applied to various output filenames (default: ${params.project}) 
      --save-config [file]         Save current command-line options to YAML file (default: options.yml)
      --publish-mode [mode]        Specify how nextflow places files in output directories 
                                   (default: symlink)
      --fastqc                     Output FastQC reports for pre and post filter/merge steps 
      --split                      Split input fastq files and process in parallel
                                   (not compatible with --demultiplexed-by index)
      --split-by                   Number of sequences per split fastq chunk (default: ${params.splitBy})
      --preprocess-only            Stop after running preprocessing steps (before denoising)
      --sample-map [mapfile]       (Optional) A tab-delimited file mapping sample names to sequence-read
                                   filenames. Paired-end runs include both forward and reverse reads. Example map:
                                   ---
                                   #sample  read1                        read2
                                   sample1  B1_S7_L001_R1_001.fastq.gz   B1_S7_L001_R2_001.fastq.gz
                                   sample2  B2_S8_L001_R1_001.fastq.gz   B2_S8_L001_R2_001.fastq.gz
                                   sample3  CL1_S2_L001_R1_001.fastq.gz  CL1_S2_L001_R2_001.fastq.gz
                                   sample4  CL2_S3_L001_R1_001.fastq.gz  CL2_S3_L001_R2_001.fastq.gz 
                                   ---

    For single-end sequencing runs:
      --single                     Specify single-ended sequencing run (required)
      --reads [file/glob/dir]      Location of sequencing read(s). For runs that have already been 
                                   demultiplexed by the sequencer, this may be a glob (e.g., '*.fastq').
                                   If a directory, files will be matched using the glob '<dir>/*.f*q*'

    For paired-end sequencing runs:  
      --paired                     Specify paired-end sequencing run (required)

    To specify the location of paired-end sequence files, the following methods are availble

    Resolve paired-end reads locations directly using globs:
      --reads [glob]               If --reads is a glob, attempt to resolve paired-end reads directly
                                   e.g., "--reads '/data/run1/*{R1,R2}*.fastq.gz'"
      --fwd [glob], --rev [glob]   Resolve forward and reverse reads directly using globs
                                   e.g., --fwd 'r1/*R1*.fastq' --rev 'r2/*R2*.fastq'

    To specify the location of paired-end reads with directories, rainbow_bridge will use the following patterns:
      <reads>/*{<r1>,<r2>}*.f*q*
      <reads>/{<fwd>,<rev>}/*{<r1>,<r2>}*.f*q*
      {<fwd>,<rev>}/*{<r1>,<r2>}*.f*q*

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
      --mate-separator [char]       Forward/reverse read mate separator 
                                    (default: '${params.mateSeparator}')
      --min-quality [num]           Minimum Phred score for sequence retention 
                                    (default: ${params.minQuality})
      --min-align-len [num]         Minimum sequence overlap when merging forward/reverse reads
                                    (default: ${params.minAlignLen})
      --min-len [num]               Minimum overall sequence length (default: ${params.minLen})

    BLAST settings:
      --blast [blastdb]             Location of BLAST database (path AND name).
                                    e.g., '/drive1/blast/custom_db', where the database files are named
                                    things like custom_db.ndb, custom_db.nhr, custom_db.nin, etc.
                                    To specify multiple BLAST databases, pass them as a list in a
                                    parameters file (see README for more information).
      --blast-taxdb [file]?         Specify NCBI taxdb.tar.gz file. 
                                    By default, the pipeline will use any taxdb files that are found
                                    alongside the existing BLAST databases. Pass this option without
                                    arguments to download taxdb from the NCBI servers. Otherwise, provide
                                    a path to a local copy of taxdb.tar.gz.
      --blast-taxa [taxa]           Limit BLAST query to the specified taxa. Multiple taxa can be passed
                                    as a comma-separated list.
      --blast-exclude-taxa [taxa]   Exclude specified taxa from BLAST search. Multiple taxa can be passed
                                    as a comma-separated list.
                                    NOTE: only one of --blast-taxa or --blast-exclude-taxa may be given.
      --blast-task [task]           Set blast+ task (default: "blastn")
      --max-query-results [num]     Maxmimum number of BLAST results to return per zOTU (default: 10)
      --percent-identity [num]      Minimum percent identity of matches to report (default: 95)
      --evalue [num]                Expectation value threshold for saving hits (default: 0.001)
      --qcov [num]                  Percent query coverage per hsp (default: 100)
      --blastn-<option> [arg?]      Pass <option> to blastn executable with optional argument

    General taxonomy operations:
      --standalone-taxonomy         Run LCA and/or insect in standalone mode (independent of pipeline)
      --ncbi-taxdump [file]         Local copy of the NCBI new_taxdump.zip archive (default: downloaded from NCBI server)
      --no-taxdump                  Suppress download of NCBI taxonomy dumps archive

    LCA taxonomy collapse:
      --lca                         Collapse assigned BLAST results by least common ancestor (LCA)
      --blast-file [file]           Blast result table (only for standalone LCA assignment)
      --seq-table [file]            zOTU table file (only for standalone LCA assignment)
      --lca-lineage [file]          Tabular file (TSV/CSV) matching taxonomic IDs (taxids) to taxonomic lineage 
      --lineage-priority            Give priority (over NCBI) to taxa in the custom lineage file
      --lca-qcov [num]              Minimum query coverage for LCA taxonomy assignment (default: 100)
      --lca-evalue [num]            Maximum e-value for LCA taxonomy refinement (default: 0.001)
      --lca-pid [num]               Minimum percent identity for LCA taxonomy assignment (default: 97)
      --lca-diff [num]              The difference between percent identities (when query coverage is
                                    identical) where species-level taxonomy is retained (default: 1)
      --lca-taxon-filter [regex]    Regular expression to filter taxa from BLAST results. 
                                    (default: '${params.lcaTaxonFilter}')
      --lca-case-insensitive        Ignore case when applying taxon filter regex.
      --lca-filter-max-qcov         Retain only BLAST results having highest query coverage.
      --dropped [str]               Placeholder text for dropped taxonomic levels (use "NA" for blank/NA)
                                    (default: '${params.dropped}')

    Insect taxonomy classification:
      --insect [classifier]         Perform taxonomy assignment using insect
                                    Accepted values of [classifier] are:
                                    - Filename to local .rds containing classifier model
                                    - One of the following (case-insensitive) primer names: 
                                      MiFish, Crust16S, Fish16S, 18SUni, 18SV4, p23S, mlCOIint, SCL5.8S
                                      (see https://github.com/shaunpwilkinson/insect#classifying-sequences)
      --insect-sequences [file]     FASTA file containing sequences to be classified (only for standalone insect clasification)
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
      --phyloseq                    Create phyloseq object (requires --lca)
      --metadata [file]             Comma or tab-separated sample metadata file (required)
      --taxonomy [tax]              Taxonomic classifaction scheme. May be pipeline-generated or user supplied
                                    (acceptable options: lca, insect, combined, <filename>; default: ${params.taxonomy})
      --tree                        Generate a phylogenetic tree (default: false)
      --optimize-tree               Attempt to optimize generation of the tree (may take a long time) (default: false)

    Data cleanup & rarefaction:
      --abundance-filter            Filter reads by minimum relative abundance
      --abundance-threshold [num]   Minimum relative abundance below which counts are set to zero
      --filter-minimum              Filter samples by minimum read count
      --min-reads [num]             Total read count below which samples will be filtered
      --rarefy                      Rarefy samples to minimum depth
      --rarefaction-method [method] Rarefaction method (available options: perm, phyloseq)
      --permutations [num]          For 'perm' rarefaction method, number of permutations
    
    Final output:
      --lca-table                   Produce final zOTU table merged with LCA taxonomy only
      --insect-table                Produce final zOTU table merged with insect taxonomy only

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
      --remove-ambiguous-indices    Removes reads with ambiguous indices in the header (i.e., not A,G,C,T)
                                    (only applies to previously-demultiplexed runs with indices in header)
      --sequences [file]            Skip demultiplexing step and use supplied FASTA 
                                    (must be in usearch/vsearch format)
      --primer-mismatch             Allowed number of mismatched primer bases 
                                    (default: ${params.primerMismatch})
      --no-primers                  Skip primer matching (ngsfilter/cutadapt) altogether. 
                                    Use with e.g., demultiplexed runs lacking primer sequences.

    Denoising and zOTU inference (usearch/vsearch options):  
      --denoiser [usearch/vsearch]  Sets the tool used for denoising & chimera removal (default: vsearch)
      --alpha [num]                 Sets the alpha parameter for the UNOISE3 algorithm (default: ${params.alpha})
      --min-abundance [num]         Minimum sequence abundance for zOTU determination; sequences below threshold will be discarded
                                    (default: ${params.minAbundance}) 
      --zotu-identity [num]         Fractional identity (0–1) for zOTU search (default: 0.97)
      --chimera-ref [file]          FASTA file to use in reference-based chimera detection

    LULU zOTU curation:
      --lulu                        Curate zOTUs using LULU
      --lulu-min-ratio-type [num]   LULU minimum ratio type (accepted values: 'min', 'avg', default: ${params.luluMinRatioType})
      --lulu-min-ratio [num]        LULU minimum ratio (default: ${params.luluMinRatio})
      --lulu-min-match [num]        LULU minimum threshold of sequence similarity to consider zOTUs as spurious (default: ${params.luluMinMatch})
                                    Choose higher values when using markers with lower genetic variation 
                                    and/or few expected PCR and sequencing errors. (default: ${params.luluMinMatch})
      --lulu-min-rc [num]           LULU minimum relative co-occurence rate (default: ${params.luluMinRc})

    Denoising and ASV inference (dada2 options):  
      --max-len [num]               Maximum overall sequence length (only for cutadapt/dada2, default: ${params.minLen})
      --plot-only                   Terminate pipeline after plotting quality profiles 
      --plot-qualities              Output quality profile plots
      --plot-errors                 Output error profile plots
      --plot-qualities-n [num]      Number of reads to sample from fastq when plotting quality profiles (default: ${params.plotQualitiesN})
      --dada-truncate [num]         Truncate reads after specified number of bases (applies to both directions, default: ${params.dadaTruncate})
      --dada-trunc-f [num]          Truncate forward reads after specified number of bases (default: ${params.dadaTruncF})
      --dada-trunc-r [num]          Truncate reverse reads after specified number of bases (default: ${params.dadaTruncR}) 
      --dada-trunc-q [num]          Truncate reads at first quality score below value (default: ${params.dadaTruncQ}) 
      --dada-max-n [num]            Discard sequences with N's over value (default: ${params.dadaMaxN}) 
      --dada-max-ee [num]           Discard reads with with higher than specified number of "expected errors" (both directions, default: ${params.dadaMaxEe}) 
      --dada-max-ee-f [num]         Discard forward reads with with higher than specified number of "expected errors" (default: ${params.dadaMaxEeF}) 
      --dada-max-ee-r [num]         Discard reverse reads with with higher than specified number of "expected errors" (default: ${params.dadaMaxEeR}) 
      --dada-remove-phix            Discard reads matching known phiX sequences
      --dada-trim-left [num]        Trim specified bases from beginning of reads (default: ${params.dadaTrimLeft})
      --dada-trim-right [num]       Trim specified bases from end of reads (default: ${params.dadaTrimRight})
      --dada-max-len [num]          Remove reads longer than specified length (default: ${params.dadaMaxLen})
      --dada-min-len [num]          Remove reads shorter than specified length (default: ${params.dadaMinLen})
      --dada-min-q [num]            Remove post-truncation reads with quality scores below value (default: ${params.dadaMinQ})
      --dada-chimera-method [mehod] Chimera detection method (default: ${params.dadaChimeraMethod})

    Resource allocation:
      --max-memory [mem]            Maximum memory available to nextflow processes, e.g., '8.GB' (default: ${params.maxMemory})
      --max-cpus [num]              Maximum cores available to nextflow processes default: ${params.maxCpus})
      --max-time [time]             Maximum time allocated to each pipeline process, e.g., '2.h' (default: ${params.maxTime})
      --max-retries [num]           The maxmimum number of times rainbow_bridge will attempt to re-execute a process
                                    that fails due to resource limitations (with increased resources for each iteration) 
                                    (default: 1)

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
