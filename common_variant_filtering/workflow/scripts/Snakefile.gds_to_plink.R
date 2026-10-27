# Convert filtered GDS to PLINK binary format
library(SeqArray)
library(SNPRelate)
library(tidyverse)

# Inputs and outputs
input_gds <- snakemake@input$gds
plink_prefix <- snakemake@params$plink_prefix

# Convert filtered GDS to PLINK format
seqGDS2BED(input_gds, plink_prefix, verbose = TRUE)