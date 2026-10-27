set.seed(19690329)

library(tidyverse)
library(magrittr)
library(SeqVarTools)
library(SNPRelate)
library(annotatr)
library(GenomicRanges)

gds_fn <- "/sc/arion/projects/load/users/patelt16/projects/adsp/adsp_gwas/20240118.adsp_pipeline_v8_master_pipeline/multiallelic_gwas/vcf_to_gds_multi/output/adsp.r4.wgs/adsp.r4.wgs.biallelic_decomposed_multiallelic.chrall.gds"
variant_annot_tsv <- "/sc/arion/projects/load/users/patelt16/projects/adsp/adsp_gwas/20240118.adsp_pipeline_v8_master_pipeline/multiallelic_gwas/vcf_to_gds_multi/output/adsp.r4.wgs/adsp.r4.wgs.biallelic_decomposed_multiallelic.gnomad.v4.0_annot.chrall.variant_annot.tsv"
gds_filtered <- "/sc/arion/projects/load/users/patelt16/projects/adsp/adsp_gwas/20240118.adsp_pipeline_v8_master_pipeline/multiallelic_gwas/all_gwas_variants/output/adsp.r4.wgs/cc.all/diagnosis_code/all/common_all/adsp.r4.wgs_cc.all.chrall.filtered.gds"
snp_filter <- NULL
impact_filter <- NULL
snp_maf <- 0.01
snp_type <- "common"

# gds_fn <- snakemake@input$gds
# variant_annot_tsv <- snakemake@input$variant_annot
# gds_filtered <- snakemake@output$gds
# snp_filter <- snakemake@params$snp_filter
# impact_filter <- snakemake@params$impact_filter
# snp_maf <- snakemake@params$snp_maf
# snp_type <- snakemake@params$snp_type

gds <- seqOpen(gds_fn)

seqResetFilter(gds)

message("Set proper variant id")
newfile <- tempfile()
file.copy(gds_fn, newfile)
gds <- seqOpen(newfile)
variant.id <- seqGetData(gds, "annotation/id")

## Ensure all variant.id are unique (issue in R4)
if (length(unique(variant.id)) != length(variant.id)) {
  variant.id <- make.unique(variant.id)
}

seqClose(gds)

setVariantID(newfile, variant.id)
gds <- seqOpen(newfile)

message("Filter for biallelic variants only")
seqSetFilter(gds, variant.sel = isSNV(gds, biallelic=TRUE))

message("Filter for coding variants using snpEff annotations")
filter_varanno <- function(x, pos) {
  x |>
    rename_at(vars(starts_with("ANN[0].")),
              funs(str_replace(., "ANN\\[0\\]\\.", ""))) |>
    dplyr::rename(variant.id = ID)
}

variant_annot <- read_tsv_chunked(variant_annot_tsv,
                                  DataFrameCallback$new(filter_varanno), chunk_size = 100000)

colnames(variant_annot)

## Extract and filter variants by MAF
variant_af_filter <- seqGetData(gds, c("variant.id", "annotation/info/AF")) |>
  as_tibble() |>
  dplyr::rename(variant.id = 1, AF = 2)

variant_af_filter <- variant_af_filter[1:1000, ]

## Annotate the variant_annot table
variant_annot_af <- variant_annot |>
  left_join(variant_af_filter, by = "variant.id") |>
  mutate(maf_filter = case_when(AF < snp_maf ~ "rare",
                                AF >= snp_maf ~ "common",
                                TRUE ~ NA))

variant_annot_af |> dplyr::count(maf_filter)

## Filter variants
variant_annot_filter <- variant_annot_af |>
  dplyr::filter(maf_filter == snp_type,
                EFFECT != "sequence_feature",
                BIOTYPE == snp_filter,
                IMPACT %in% impact_filter) |>
  mutate(START = POS, END = POS + 1) |>
  pull(variant.id)

nrow(variant_annot_filter)

seqSetFilter(gds, variant.id = unique(variant_annot_filter))

## Filter for exonic variants
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# txs <- exonsBy(txdb, "gene")

# seqlevels(txs) <- gsub("chr", "", seqlevels(txs))

# seqSetFilter(gds, txs)

message("Saving filtered GDS file")

seqExport(gds, gds_filtered, verbose = TRUE)

seqClose(gds)

sessionInfo()

# showfile.gds(closeall=TRUE)
