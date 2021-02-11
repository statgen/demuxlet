# demuxlet
Genetic multiplexing of barcoded single cell RNA-seq

### Citation

Demuxlet has been published: https://www.nature.com/articles/nbt.4042

If you find it useful, please cite: Kang et al., Nature Biotechnology 2017.

### Tips for running

* Set `--alpha 0 --alpha 0.5`, which assumes the expected proportion of 50% genetic mixture from two individuals, to get better estimates of doublets.
* Set `--group-list` to a list of barcodes (i.e. barcodes.tsv from 10X) to speed things up and only get demultiplexing for cells called by other methods
* To reproduce the results presented in Figure 2 of the demuxlet paper, please go to: https://github.com/yelabucsf/demuxlet_paper_code/tree/master/fig2 to download the vcf and the outputs of demuxlet.

### Introduction

**_demuxlet_** is a software tool to deconvolute sample identity and identify multiplets when multiple samples are pooled by barcoded single cell sequencing.
**_demuxlet_** takes (1) a SAM/BAM/CRAM file produced by the standard 10x sequencing platform, or any other barcoded single cell RNA-seq (with proper --tag-UMI and --tag-group) options (2) a VCF/BCF file containing the genotype (GT), posterior probability (GP), or genotype likelihood (GL) to assign each barcode to a specific sample (or a pair of samples) in the VCF file. 

### Installing demuxlet

Before installing demuxlet, you need to install [htslib](https://github.com/samtools/htslib) in the same directory you want to install demuxlet (i.e. demuxlet and htslib should be siblings).  
**NOTE** htslib 1.11 is not supported for now - use earlier releases (e.g. 1.10.x)

After installing htslib, you can clone the current snapshot of this repository to install as well

<pre>
$ git clone https://github.com/statgen/demuxlet.git
$ cd demuxlet
$ autoreconf -vfi
$ ./configure  (with additional options such as --prefix)
$ make
$ make install (may require root privilege)
</pre>

### Using demuxlet

demuxlet uses a self-documentation utility. You can run each utility with -man or -help option to see the command line usages.

<pre>
$ ./demuxlet          (for short usage)
$ ./demuxlet -help    (for detailed usage)
</pre>

The detailed usage is also pasted below.

<pre>
Options for input SAM/BAM/CRAM
  --sam           [STR: ]             : Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed
  --tag-group     [STR: CB]           : Tag representing readgroup or cell barcodes, in the case to partition the BAM file into multiple groups. For 10x genomics, use CB
  --tag-UMI       [STR: UB]           : Tag representing UMIs. For 10x genomiucs, use UB

Options for input VCF/BCF
  --vcf           [STR: ]             : Input VCF/BCF file, containing the individual genotypes (GT), posterior probability (GP), or genotype likelihood (PL)
  --field         [STR: GP]           : FORMAT field to extract the genotype, likelihood, or posterior from
  --geno-error    [FLT: 0.01]         : Genotype error rate (must be used with --field GT)
  --min-mac       [INT: 1]            : Minimum minor allele frequency
  --min-callrate  [FLT: 0.50]         : Minimum call rate
  --sm            [V_STR: ]           : List of sample IDs to compare to (default: use all)
  --sm-list       [STR: ]             : File containing the list of sample IDs to compare

Output Options
  --out           [STR: ]             : Output file prefix
  --alpha         [V_FLT: ]           : Grid of alpha to search for (default is 0.1, 0.2, 0.3, 0.4, 0.5)
  --write-pair    [FLG: OFF]          : Writing the (HUGE) pair file
  --doublet-prior [FLT: 0.50]         : Prior of doublet
  --sam-verbose   [INT: 1000000]      : Verbose message frequency for SAM/BAM/CRAM
  --vcf-verbose   [INT: 10000]        : Verbose message frequency for VCF/BCF

Read filtering Options
  --cap-BQ        [INT: 40]           : Maximum base quality (higher BQ will be capped)
  --min-BQ        [INT: 13]           : Minimum base quality to consider (lower BQ will be skipped)
  --min-MQ        [INT: 20]           : Minimum mapping quality to consider (lower MQ will be ignored)
  --min-TD        [INT: 0]            : Minimum distance to the tail (lower will be ignored)
  --excl-flag     [INT: 3844]         : SAM/BAM FLAGs to be excluded

Cell/droplet filtering options
  --group-list    [STR: ]             : List of tag readgroup/cell barcode to consider in this run. All other barcodes will be ignored. This is useful for parallelized run
  --min-total     [INT: 0]            : Minimum number of total reads for a droplet/cell to be considered
  --min-uniq      [INT: 0]            : Minimum number of unique reads (determined by UMI/SNP pair) for a droplet/cell to be considered
  --min-snp       [INT: 0]            : Minimum number of SNPs with coverage for a droplet/cell to be considered
</pre>

### Interpretation of output files

**_demuxlet_** generates multiple output file, such as `[prefix].best`, `[prefix].sing`, `[prefix].sing2`, and optionally `[prefix].pair` (with `--write-pair` argument). Each file contains the following information
* The `[prefix].best` file contains the best guess of the sample identity, with detailed statistics to reach to the best guess
* The `[prefix].sing` file contains the statistics for matching each cell with each possible sample.
* The `[prefix].sing2` file contains the statistics similar information to the previous one, but generated for sanity checking of the `[prefix].pair` results.
* The `[prefix].pair` file contains the statistics for matching each cell with each possible configuration of doublet. 

The `[prefix].best` file contains the following 22 columns.
 1. BARCODE - Cell barcode for the cell that is being assigned in this row
 2. RD.TOTL - The total number of reads overlapping with variant sites for each droplet.
 3. RD.PASS - The total number of reads that passed the quality threshold, such as mapping quality, base quality. 
 4. RD.UNIQ - The total number of UMIs that passed the quality threshold. If a UMI is observed in a single variant multiple times, it won't be counted more. If a UMI is observed across multiple variants, it will be counted as different.
 5. N.SNP   - The total number of variants overlapping with any read in the droplet.
 6. BEST    - The best assignment for sample ID.
    * For singlets, SNG-<sample ID>
    * For doublets, DBL-<sample ID1>-<sampleID2>-<mixture rate>
    * For ambiguous droplets, , AMB-<best-singlet-sampleID>-<next-best-singlet-sampleID>-<doublet ID1/ID2>)
 7. SNG.1ST - The best singlet assignment for sample ID
 8. SNG.LLK1 - The log(likelihood that the ID from SNG.1ST is the correct assignment)       
 9. SNG.2ND - The next best singlet assignment for sample ID
 10. SNG.LLK2 - The log(likelihood that the ID from SNG.2ND is the correct assignment)        
 11. SNG.LLK0 - The log-likelihood from allele frequencies only      
 12. DBL.1ST - The sample ID that is most likely included if the assignment is a doublet
 13. DBL.2ND - The sample ID that is next most likely included ifthe assignment is a doublet
 14. ALPHA   - % Mixture Proportion
 15. LLK12   - The log(likelihood that the ID is a doublet)
 16. LLK1    - The log(likelihood that the ID from DBL.1ST is the correct singlet assignment)
 17. LLK2    - The log(likelihood that the ID from DBL.2ND is the correct singlet assignment)
 18. LLK10   - The log(likelihood that the ID from DBL.1ST is one of the doublet, and the other doublet identity is calculated from allele frequencies only)   
 19. LLK20   - The log(likelihood that the ID from DBL.2ND is one of the doublet, and the other doublet identity is calculated from allele frequencies only)   
 20. LLK00   - The log(likelihood that the droplet is doublet, but both identities are calculated from allele frequencies only)
 21. PRB.DBL - Posterior probability of the doublet assignment
 22. PRB.SNG1 - Posterior probability of the singlet assignment when excluding all possible doublets
    
