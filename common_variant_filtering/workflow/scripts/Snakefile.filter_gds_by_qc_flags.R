# Filter GDS based on QC flags from variant_filterlist
library(SeqArray)
library(tidyverse)
library(magrittr)

# Inputs and outputs
# filterlist_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/common_all/report/adsp.r4.wgs_cc.all.chrall.variant_qc_flags_filterlist.tsv"
# output_gds <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/common_all/adsp.r4.wgs_cc.all.chrall.filtered_qc.gds"

input_gds <- snakemake@input$gds
filterlist_tsv <- snakemake@input$variant_filterlist
output_gds <- snakemake@output$filtered_gds

# Read the variant filterlist
cat("Reading variant filterlist from:", filterlist_tsv)
variant_flags <- read_tsv(filterlist_tsv, col_types=cols(.default="c"))

cat("Variant flags:")
print(dplyr::count(variant_flags, FILTER_FLAG))

# Open input GDS
gds <- seqOpen(input_gds)

# Numeric IDs and CPRA IDs in the same GDS variant order
gds_ids  <- seqGetData(gds, "variant.id")
gds_cpra <- seqGetData(gds, "annotation/id")
# gds_cpra <- seqGetData(gds, "$chrom_pos_allele")

# PASS CPRA strings from your flags table
pass_cpra <- variant_flags |>
  dplyr::filter(FILTER_FLAG == "PASS") |>
  dplyr::pull(variant.id)

# Map CPRA matches back to internal numeric IDs
variants_to_keep <- gds_ids[gds_cpra %in% pass_cpra]

message("Keeping ", length(variants_to_keep), " variants")
stopifnot(length(variants_to_keep) > 0L)

# Filter to keep only high quality variants
seqSetFilter(gds, variant.id = variants_to_keep)

message("Save QC filtered GDS")
seqExport(gds, output_gds, verbose = TRUE)

# Close gds
seqClose(gds)
