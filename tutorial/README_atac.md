# Using demuxlet for single cell ATACseq
## The following was tested and applies to the 10x Genomics scATAC platform 
Demuxlet works for single cell ATACseq data out of the box. In most cases you can use the same `.vcf` you would use
for the RNAseq. However, we also recommend using a broader `.vcf` file which contains regions from the entire genome. 
Try to stick to the best practices of `.vcf` filtering as described in the `README_vcf.md`. 

 
There are slight adjustments that have to be made with the arguments for the running command.  
Specifically, there is no need to specify `--tag-UMI` because 10x platform does not have UMI sequence associated with any read.  
