# Merge chromosome VCF files into a single GDS file
library(SeqArray)
library(tidyverse)

# Get input and output files
chr_vcf_files <- snakemake@input$chr_vcf
output_gds <- snakemake@output$merged_gds

cat("Merging", length(chr_vcf_files), "VCF files into GDS...")

## Merge vcfs and convert to gds using LZ4 compression instead of LZMA
seqVCF2GDS(unlist(chr_vcf_files), output_gds, parallel = 80, storage.option = "LZ4_RA", verbose = TRUE)

cat("GDS file created successfully:", output_gds)