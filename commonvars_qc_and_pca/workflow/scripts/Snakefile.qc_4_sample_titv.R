set.seed(19690329)

## Load packages
library(tidyverse)
library(magrittr)
library(SeqVarTools)
library(SNPRelate)
library(readr)

## Inputs
# gds_fn <- "/sc/arion/projects/load/users/patelt16/projects/adsp/adsp_gwas/20240118.adsp_pipeline_v8_master_pipeline/multiallelic_gwas/vcf_to_gds_multi/output/adsp.r4.wgs/adsp.r4.wgs.biallelic_decomposed_multiallelic.chrall.gds"
# sample_pheno_tsv <- "/sc/arion/projects/load/users/patelt16/projects/adsp/adsp_gwas/20240118.adsp_pipeline_v8_master_pipeline/phen/20250519.multiallelic_run/20250519.ADSP_Umbrella_ng00067.v10.wgs.cc.all.28063.tsv"
# variant_missingness_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/adsp.r4.wgs_cc.all.variant_metrics_missingness.tsv"
# sample_titv_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/adsp.r4.wgs_cc.all.sample_metrics.tsv"
# qc_summary_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/qc_metrics.tsv"
# qc_report_pdf <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/qc_report.pdf"
# snp_maf <- 0.01
# missing.rate <- 0.1
# refp.super_pops <- "all"

gds_fn <- snakemake@input$gds
sample_pheno_tsv <- snakemake@input$sample_pheno_tsv
variant_missingness_tsv <- snakemake@input$variant_missingness_tsv
sample_titv_tsv <- snakemake@output$sample_titv_tsv
# qc_summary_tsv <- snakemake@output$qc_summary_tsv
qc_report_pdf <- snakemake@output$qc_report_pdf
snp_maf <- snakemake@params$snp_maf
missing.rate <- snakemake@params$missing_rate
refp.super_pops <- snakemake@params$refp_super_pops

## Script
`%notin%` <- Negate(`%in%`)

gds <- seqOpen(gds_fn)
# seqSetFilterChrom(gds, include = 22, verbose=TRUE)

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

## Filter for variant missing rate
variant_metrics_filter <- read_tsv(variant_missingness_tsv) |>
  dplyr::filter(!is.na(pass_missingness))

## Filter for samples in pheno file and biallelic variants
message("Filter for samples in pheno file and biallelic variants only")
seqSetFilter(gds, sample.id = sample_pheno$sample.id, variant.sel = isSNV(gds, biallelic=TRUE))

seqSetFilter(gds, variant.id = variant_metrics_filter$variant.id)

## Calculate sample Ti/Tv ratio
message("Ti/Tv ratio should be approximately 2.1 for whole genome sequencing and 2.8 for whole exome sequencing!")
sample.titv <- titv(gds, by.sample=TRUE)

sample.metrics <- tibble(
  sample.id = seqGetData(gds, "sample.id"),
  titv_ratio = sample.titv,
  zscore = as.numeric(scale(sample.titv)),
  pass_titv = abs(scale(sample.titv)[,1]) <= 3
)

write_tsv(sample.metrics, sample_titv_tsv)

## Plot sample hethom rate
pdf(qc_report_pdf, width = 10, height = 7.5)

ggplot(sample.metrics, aes(x = titv_ratio)) +
  geom_histogram(bins = 60, fill = "#3182bd", color = "white") +
  geom_vline(xintercept = c(2.1, 2.8), color = "red", linetype = "dashed") +
  labs(x="Ti/Tv Ratio", y="Samples",
       title="Sample Transition/Transversion Ratio") +
  theme_minimal(base_size = 14)

dev.off()

seqClose(gds)
# showfile.gds(closeall=TRUE)

sessionInfo()

