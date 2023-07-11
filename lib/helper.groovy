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
    Usage: eDNAFlow.nf [options]

    General options:
      --prefix [prefix]          Project prefix, applied to output filenames as "sample ID" 
                                 (default: 'seq')
      --barcode [file]           (required) Barcode file. Format must match OBITools requirements
                                 (see https://pythonhosted.org/OBITools/scripts/ngsfilter.html)
                                 To denote multiple barcode files, this may be a glob, but it must
                                 be surrounded by quotations (e.g. 'barcode*.tab'). 
                                 For previously-demultiplexed datasets, use ':' for barcode pairs
                                 Primer sequences are still used for primer-mismatch comparisons
      --publish-mode [mode]      Specify how nextflow places files in output directories 
                                 (default: symlink)
      --fastqc                   Output FastQC reports for pre and post filter/merge steps 

    For single-end sequencing runs:
      --single                   Specify single-ended sequencing run (required)
      --reads [file/glob]        Location of sequencing read(s). For runs that have already been 
                                 demultiplexed by the sequencer, this may be a glob (e.g., '*.fastq')

    For paired-end sequencing runs:  
      --paired                   Specify paired-end sequencing run (required)

    To specify the location of paired-end sequence files, the following file patterns are followed:
      <base-dir>/*<r1>|<r2>*.fastq
      <base-dir>/<fwd>|<rev>/*<r1>|<r2>*.fastq

    Options:
      --base-dir [dir]           Base directory of sequence read subdirectories (default: .) 
      --fwd [dir]                (Optional) forward reads directory (default: ${params.fwd})
      --rev [dir]                (Optional) reverse reads directory (default: ${params.rev})
      --r1 [pattern]             Pattern distinguishing forward read files (default: ${params.r1})
      --r2 [pattern]             Pattern distinguishing reverse read files (default: ${params.r2})

    Length, quality, and merge settings:
      --min-quality [num]        Minimum Phred score for sequence retention 
                                 (default: ${params.minQuality})
      --min-align-len [num]      Minimum sequence overlap when merging forward/reverse reads
                                 (default: ${params.minAlignLen}
      --min-len [num]            Minimum overall sequence length (default: ${params.minLen})

    BLAST settings (one or more of --blast-db or --custom-db is required):
      --blast-db [dir]           Location of local BLAST nucleotide (nt) database *directory*
                                 (do NOT include the final "/nt")
                                 May be specified using the \$BLASTDB environment variable
                                 (current value: ${get_env("BLASTDB")})
      --custom-db [dir]          (optional) Path to custom BLAST database *directory*
                                 Passing --custom-db and --blast-db values that point to directories
                                 with the same name (e.g., /dir1/blast and /dir2/blast) will result 
                                 in an error!
      --custom-db-name [name]    Name of custom BLAST database (i.e., basename of .ndb, etc. files)
      --blast-task [task]        Set blast+ task (default: "blastn")
      --max-query-results [num]  Maxmimum number of BLAST results to return per zOTU (default: 10)
      --percent-identity [num]   Minimum percent identity of matches to report (default: 95)
      --evalue [num]             Expectation value threshold for saving hits (default: 0.001)
      --qcov [num]               Percent query coverage per hsp (default: 100)

    Taxonomy assignment:
      --assign-taxonomy          Perform final taxonomy assignment & LCA collapse 
      --lca-qcov [num]           Minimum query coverage for LCA taxonomy assignment (default: 100)
      --lca-pid [num]            Minimum percent identity for LCA taxonomy assignment (default: 97)
      --lca-diff [num]           The difference between percent identities (when query coverage is
                                 identical) where species-level taxonomy is retained (default: 1)
      --insect [classifier]      Perform taxonomy assignment using insect
                                 Accepted values of [classifier] are:
                                 - Filename to local .rds containing classifier model
                                 - One of the following (case-insensitive) primer names: 
                                   MiFish, Crust16S, Fish16S, 18SUni, 18SV4, p23S, mlCOIint, SCL5.8S
                                   (see https://github.com/shaunpwilkinson/insect#classifying-sequences)
                                 

    Demultiplexing and sequence matching:
      --illumina-demultiplexed   Sequencing run has already been demultiplexed
      --remove-ambiguous-tags    Removes reads with ambiguous tags in the header (i.e., not A,G,C,T)
                                 (only applies to previously-demultiplexed runs with tags in header)
      --demuxed-fasta [file]     Skip demultiplexing step and use supplied FASTA 
                                 (must be in usearch/vsearch format)
      --demuxed-example          Spit out example usearch/vsearch demultiplexed FASTA format
      --demux-only               Stop after demultiplexing and splitting raw reads
      --primer-mismatch          Allowed number of mismatched primer bases 
                                 (default: ${params.primerMismatch})

    Denoising and zOTU inference:
      --min-abundance [num]      Minimum zOTU abundance; zOTUs below threshold will be discarded
                                 (default: ${params.minAbundance}) 
      --denoiser [tool/path]     Sets the tool used for denoising & chimera removal
                                 accepted options: usearch/usearch32 (equivalent), vsearch, 
                                 path to usearch64 executable (default: usearch)

    LULU zOTU curation:
      --lulu-min [num]           Minimum threshold of sequence similarity to consider zOTUs as spurious.
                                 Choose higher values when using markers with lower genetic variation 
                                 and/or few expected PCR and sequencing errors. (default: ${params.luluMin})

    Resource allocation:
      --max-memory [mem]         Maximum memory available to nextflow processes, e.g., '8.GB' (default: ${params.maxMemory})
      --max-cpus [num]           Maximum cores available to nextflow processes default: ${params.maxCpus})
      --max-time [time]          Maximum time allocated to each pipeline process, e.g., '2.h' (default: ${params.maxTime})

    Singularity options:
      --bind-dir [dir]           Space-separated list of directories to bind within singularity images
                                 (must be surrounded by quotations if more than one directory)
                                 Note: singularity will attempt to auto-bind all provided host paths
                                 so this option may not be necessary, but try it if you're getting 
                                 weird "file not found" types of errors
      --singularity-cache [dir]  Location to store singularity images. May also be specified
                                 with the environment variable \$NXF_SINGULARITY_CACHEDIR.
                                 (current value: ${get_env("NXF_SINGULARITY_CACHEDIR")})
		""".stripIndent())
  }
}
