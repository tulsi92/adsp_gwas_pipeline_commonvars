## Load packages
library(readr)
library(tidyverse)
library(dplyr)
library(magrittr)

## Inputs
# variant_annot_tsv <- "output/adsp.r4.wgs/adsp.r4.wgs.decomposed_multiallelic_variant_quality_metrics.chr21.tsv"
# variant_filterlist_tsv <- "output/adsp.r4.wgs/adsp.r4.wgs.variant_filterlist.chr21.tsv"
# variant_annot_abhet_allele_tsv <- "output/adsp.r4.wgs/adsp.r4.wgs.variant_annot_multiallelic_decomp.chr21.tsv"

variant_annot_tsv <- snakemake@input$variant_annot_tsv
variant_annot_abhet_allele_tsv <- snakemake@output$variant_annot_abhet_allele_tsv
variant_filterlist_tsv <- snakemake@output$variant_filterlist_tsv

## Read in variant quality flags
variant_annot <- read.delim(variant_annot_tsv)

variant_annot_filt <- variant_annot |>
  dplyr::filter(VFLAGS_One_subgroup == 0) |>
  mutate(variant.id = paste(CHROM, POS, A1, A2, sep = ":"))

variant_annot_filt |>
  dplyr::count(ABHet_One_subgroup)

variant_annot_abhet <- variant_annot_filt |> 
  dplyr::select(variant.id, CHROM, POS, VFLAGS_One_subgroup, ABHet_One_subgroup)

variant_annot_abhet_allele <- variant_annot_abhet |>
  ## Split ABHet_One_subgroup into separate values
  mutate(ABHet_One_subgroup_split = strsplit(as.character(ABHet_One_subgroup), ",")) |>
  ## Group by CHROM/POS to ensure proper indexing for duplicates
  group_by(CHROM, POS) |>
  ## Create new column for index based on row number within each group
  mutate(ABHet_One_subgroup_allele = map2_chr(ABHet_One_subgroup_split, row_number(), ~ .x[.y])) |>
  ## Remove spaces after commas in ABHet_One_subgroup
  mutate(ABHet_One_subgroup = gsub(" ", "", ABHet_One_subgroup)) |>
  ## Remove temporary split column
  dplyr::select(variant.id, CHROM, POS, VFLAGS_One_subgroup, ABHet_One_subgroup, ABHet_One_subgroup_allele) |>
  ungroup()

head(variant_annot_abhet_allele)

variant_annot_abhet_allele |>
  write_tsv(variant_annot_abhet_allele_tsv)

variant_filterlist <- variant_annot_abhet_allele |>
  dplyr::filter(ABHet_One_subgroup_allele > 0.25 & ABHet_One_subgroup_allele < 0.75) |>
  dplyr::select(CHROM, POS) |>
  write_tsv(variant_filterlist_tsv, col_names = FALSE)

