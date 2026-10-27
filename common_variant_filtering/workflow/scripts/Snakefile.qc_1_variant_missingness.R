set.seed(19690329)

## Load packages
library(tidyverse)
library(magrittr)
library(SeqVarTools)
library(SNPRelate)
library(readr)

## Inputs
# gds_fn <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/test/adsp.r4.wgs_cc.all.chrall.filtered.gds"
# sample_pheno_tsv <- "/sc/arion/projects/load/users/patelt16/projects/adsp/adsp_gwas/20240118.adsp_pipeline_v8_master_pipeline/phen/20250519.multiallelic_run/20250519.ADSP_Umbrella_ng00067.v10.wgs.cc.all.28063.tsv"
# variant_missingness_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/adsp.r4.wgs_cc.all.variant_metrics_missingness.tsv"
# qc_summary_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/qc_metrics.tsv"
# qc_report_pdf <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/qc_report.pdf"
# snp_maf <- 0.01
# missing.rate <- 0.1
# refp.super_pops <- "all"

gds_fn <- snakemake@input$gds
sample_pheno_tsv <- snakemake@input$sample_pheno_tsv
variant_missingness_tsv <- snakemake@output$variant_missingness_tsv
# qc_summary_tsv <- snakemake@output$qc_summary_tsv
qc_report_pdf <- snakemake@output$qc_report_pdf
snp_maf <- snakemake@params$snp_maf
missing.rate <- snakemake@params$missing_rate
refp.super_pops <- snakemake@params$refp_super_pops

## Script
`%notin%` <- Negate(`%in%`)

gds <- seqOpen(gds_fn)

## Sample pheno file
sample_pheno <- read_tsv(sample_pheno_tsv)

if (length(refp.super_pops) > 0 & 'super_pop' %notin% colnames(sample_pheno)) {
  sample_pheno <- sample_pheno |> mutate(super_pop = case_when(population == "NHW" ~ "EUR",
                                              population == "AA" ~ "AFR",
                                              population == "HI" ~ "AMR",
                                              is.na(population) ~ NA,
                                              TRUE ~ "OTHER"))
}

if (length(refp.super_pops) > 0 & refp.super_pops != "all" & 'super_pop' %in% colnames(sample_pheno)) {
  refp.super_pops_all <- unlist(strsplit(refp.super_pops, split = "_"))
  sample_pheno <- sample_pheno |>
    dplyr::filter(super_pop %in% c(refp.super_pops_all))
}

seqResetFilter(gds)

## Filter for samples in pheno file and biallelic variants
message("Filter for samples in pheno file and biallelic variants only")
seqSetFilter(gds, sample.id = sample_pheno$sample.id, variant.sel = isSNV(gds, biallelic=TRUE))

## Calculate variant missing rate
message("Variant-level missing genotype rate filtering")
variant.missing.rate <- missingGenotypeRate(gds, margin = "by.variant")

variant.metrics <- tibble(
  variant.id = seqGetData(gds, "annotation/id"),
  missing_rate = variant.missing.rate,
  pass_missingness = variant.missing.rate <= missing.rate
)

write_tsv(variant.metrics, variant_missingness_tsv)

## Plot variant missing rate
pdf(qc_report_pdf, width = 10, height = 7.5)

ggplot(variant.metrics, aes(x = missing_rate)) +
  geom_histogram(bins = 100, fill = "#2c7fb8", color = "white") +
  geom_vline(xintercept = missing.rate, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Variant-level Missing Genotype Rate",
       x = "Missing Genotype Rate",
       y = "Count") +
  theme_bw(base_size = 14)

dev.off()

seqClose(gds)
# showfile.gds(closeall=TRUE)

sessionInfo()
