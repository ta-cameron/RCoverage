# RCoverage
This package allows you to generate read coverage maps of paired-end Illumina sequencing data, including the entire length of the inferred fragment. Most software (eg, IGV) is designed for use with Eukaryotic genomes and will only display coverage of the actual reads... which makes perfect sense when you have the ambiguity of introns to worry about. However, for bacterial RNA-seq data, introns are not a concern and mapping the entire read + inferred fragment region gives a better idea of the actual RNA fragments present in the cell. This can be particularly important when examining the effects of RNA processing genes and sRNAs.

# Instructions
1. Load the functions in the RCoverage.R file.
2. See example plots.R for usage examples.
