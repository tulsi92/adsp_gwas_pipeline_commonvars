# install.packages("GenomeInfoDbData", repos="https://cran.rstudio.com/")
library(SeqArray)
library(GenomeInfoDbData)

input.vcf <- snakemake@input
output.gds <- as.character(snakemake@output)

## Merge vcfs and convert to gds using LZ4 compression instead of LZMA
seqVCF2GDS(unlist(input.vcf), output.gds, parallel = 80, storage.option = "LZ4_RA", verbose = TRUE)

# seqSummary(output.gds)
