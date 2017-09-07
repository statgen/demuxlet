# Tutorial for demuxlet ##

## 1. Introduction ##

### 1.1. Overview ###

This tutorial provides streamlined instructions for using the tool `demuxlet`. For a more detailed description of all of the options available to use with `demuxlet`, please refer to the README. `demuxlet` has several dependencies and this tutorial has been updated to use a docker instance to minimize compatibility issues. For detailed compilation instructions, see the `demuxlet` README.

`demuxlet` is a software tool to deconvolute sample identity and identify multiplets when multiple samples are pooled by barcoded single cell sequencing. `demuxlet` requires the following input files:

1. a SAM/BAM/CRAM file produced by the standard 10x sequencing platform, or any other barcoded single cell RNA-seq (with proper `--tag-UMI` and `--tag-group`) options
2. a VCF/BCF file containing the genotype (GT), posterior probability (GP), or genotype likelihood (GL) to assign each barcode to a specific sample (or a pair of samples) in the VCF file.

### 1.2. Additional resources ###

The README for demuxlet is available [here](https://github.com/statgen/demuxlet)

If you have questions about using demuxlet or suggestions for future releases, please contact jimmie.ye@ucsf.edu.

## 2. Getting Started ##

### 2.1. Installing docker ###

First, get docker for whatever platform you feel comfortable with: https://www.docker.com/community-edition#/download.

### 2.2. Running demuxlet from docker ###

We have created a docker container at `yimmieg/demuxlet` to run `demuxlet` through docker.

You can run it with 

```
$ docker run -v path/to/tutorial:/data yimmieg/demuxlet --sam /data/$bam --vcf /data/$vcf --field $(GT or GP or PL) --out /data/$filename
```

Add bam file name for `$bam` and vcf file name for `$vcf`. Use `<(zcat $vcf)` or `<(gzcat $vcf)` if vcf file is compressed

The options for `--field` are individual genotypes (GT), posterior probability (GP), or genotype likelihood (PL). If using GT option for `--field`, you must include `--geno-error`, which is the genotype error rate.

### 2.3. demuxlet output ###

The demuxlet software produces three output files.

1. `[prefix].best`
The .best file contains the assignments of the best sample identity (singlet: `SNG-<sample name>`; doublet: `DBL-<sample IDs>`; ambiguous: `AMB-< >`) in the `BEST` column for each cell barcode identified in the `BARCODE` column along with details of the statistics used to determine the best identity.
2. The `[prefix].single` file contains the statistics for matching each cell with each possible sample

For complete descriptions of the columns in each output file, please see the demuxlet README.

## 3. Analyzing the sample dataset ##

### 3.1 Create directory and download datasets ###

Now, let's first download the data we need. We are now providing a one-stop-shop to download all of the data you need for the tutorial: https://ucsf.box.com/s/vg1bycvsjgyg63gkqsputprq5rxzjl6k. After downloading and unzipping (it's BIG), you should have a directory called `demuxlet.tutorial`.

### 3.2 Run demuxlet ###

```
$ cd demuxlet.tutorial
$ docker run -v ./:/data yimmieg/demuxlet --sam /data/jurkat_293t_downsampled_n500_full_bam.bam --vcf /data/jurkat_293t_exons_only.vcf --field GT --out data/jurkat_293t_demuxlet
```

### 3.3 Compare the called genotypes vs transcriptome data ###

In this analysis, we will use R to produce a t-SNE (t-Distributed Stochastic Neighbor Embedding) plot of the cells from the 293T:Jurkat 10x experiment with the cells colored by the assignments from the demuxlet pipeline.

```
library(ggplot2);
library(data.table);

## let's read in the barcodes
tsne <- fread("analysis_csv/tsne/projection.csv");
demuxlet <- fread("jurkat_293t_demuxlet.best");

## let's filter for the barcodes that we sampled

df <- data.frame(tsne1=tsne$"TSNE-1"[na.omit(match(demuxlet$BARCODE,tsne$Barcode))], tsne2=tsne$"TSNE-2"[na.omit(match(demuxlet$BARCODE,tsne$Barcode))], doublet=sapply(demuxlet$BEST,function(x){strsplit(x,"-")[[1]][[1]]}),
cell.type=sapply(demuxlet$BEST,function(x){strsplit(x,"-")[[1]][[2]]}))

ggplot(aes(tsne1,tsne2,color=cell.type),data=df)+geom_point()
```


## 4. Source of 293T and Jurkat VCF file. ##
1. 293T VCF
Source website:
http://hek293genome.org/v2/data.php
Source file:
http://bioinformatics.psb.ugent.be/downloads/genomeview/hek293/SNP/293T_RTG.vcf.gz
2. Jurkat VCF
Source website: https://zenodo.org/record/400615#.WYIh7IQrLIV
Source file: https://zenodo.org/record/400615/files/jurkat_final_variant_calls.tar.gz
3. 293T:jurkat VCF file generation
We used the CrossMap tool to liftover the 293T vcf file from hg18 to hg19. The tetraploid genotype for the Jurkat vcf was collapsed to a diploid genotype before being merged with the 293T vcf file and the resulting file was filtered to contain only the exon positions.
    
