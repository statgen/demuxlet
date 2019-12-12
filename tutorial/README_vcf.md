# Filtering of .vcf
These are some recommendations how to filter your `.vcf` file for demuxlet and freemuxlet tools.  


1. Adhere to Broad Institute’s GATK best practices whenever possible.  
2. Apply all or some of the additional filters:  
    Classified as common in the 1000 genomes project.   
    Called with high variant confidence (QD > 10.0)  
    Called in every one of your samples.
    Minor Allele Frequency (MAF 0.01)
    Minor Allele Counts, if applicable (MAC 1)  
3. **For RNA:** Filter for exonic variants.  
    **For ATAC** Filter out blacklisted variants (e.g. https://www.encodeproject.org/annotations/ENCSR636HFF/)
    
## Examples  
We attach some `.bed` files with positions of SNPs which have been filtered in our lab and confirmed to perform well. They are based on hg19.  
You can use these sites to filter your `.vcf` as an easy initial step. Please note that you may have a very low overlap between these and your variants,  
depending on the genotyping and variant calling strategies. 