set.seed(19620321)

library(tidyverse)
library(magrittr)
library(SeqVarTools)
library(SNPRelate)
library(readr)

# gds_fn <- "/sc/arion/projects/load/users/patelt16/projects/adsp/adsp_gwas/20240118.adsp_pipeline_v8_master_pipeline/multiallelic_gwas/vcf_to_gds_multi/output/adsp.r4.wgs/adsp.r4.wgs.biallelic_decomposed_multiallelic.chrall_pca_maf0.01.gds"
# sample_pheno_tsv <- "/sc/arion/projects/load/users/patelt16/projects/adsp/adsp_gwas/20240118.adsp_pipeline_v8_master_pipeline/phen/20250519.multiallelic_run/20250519.ADSP_Umbrella_ng00067.v10.wgs.cc.all.28063.tsv"
# variant_qc_summary_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/common_all/report/adsp.r4.wgs_cc.all.chrall.variant_qc_summary.tsv"
# pca_snp <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/common_all/adsp.r4.wgs_cc.all.chrall_pca.snp"
# pca_smp <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/common_all/adsp.r4.wgs_cc.all.chrall_pca.smp"
# pca_maf <- 0.05
# snp_missing_rate <- 0.05
# refp.super_pops <- "all"

gds_fn <- snakemake@input$gds
sample_pheno_tsv <- snakemake@input$sample_pheno_tsv
variant_qc_summary_tsv <- snakemake@input$variant_qc_summary_tsv
pca_snp <- snakemake@output$pca_snp
pca_smp <- snakemake@output$pca_smp
pca_qc_summary_tsv <- snakemake@output$pca_qc_summary_tsv
pca_maf <- snakemake@params$pca_maf
snp_missing_rate <- snakemake@params$missing_rate
refp.super_pops <- snakemake@params$refp_super_pops

`%notin%` <- Negate(`%in%`)

gds <- seqOpen(gds_fn)

sample_pheno <- read_tsv(sample_pheno_tsv)

if (length(refp.super_pops) > 0 & 'super_pop' %notin% colnames(sample_pheno)) {
  sample_pheno <- sample_pheno %>% mutate(super_pop = case_when(population == "NHW" ~ "EUR",
                                             population == "AA" ~ "AFR",
                                             population == "HI" ~ "AMR",
                                             is.na(population) ~ NA,
                                             TRUE ~ "OTHER"))
}

if (length(refp.super_pops) > 0 & refp.super_pops != "all" & 'super_pop' %in% colnames(sample_pheno)) {
  refp.super_pops_all <- unlist(strsplit(refp.super_pops, split = "_"))
  sample_pheno %<>% dplyr::filter(super_pop %in% c(refp.super_pops_all))
}

## Filter for sample and variant quality metrics
seqResetFilter(gds)

df.pca.metrics <- data.frame(Threshold = numeric(0), PCA_N_Variants = numeric(0), PCA_N_Samples = numeric(0))
df.pca.metrics["Total", ] <- c(Threshold = "null", PCA_N_Variants = length(seqGetData(gds, "variant.id")), PCA_N_Samples = length(seqGetData(gds, "sample.id")))

seqSetFilter(gds, sample.id = sample_pheno$sample.id, variant.sel = isSNV(gds, biallelic=TRUE))

df.pca.metrics["Biallelic_and_analysis_samples", ] <- c("null", length(seqGetData(gds, "variant.id")), length(seqGetData(gds, "sample.id")))

# seqSetFilterCond(gds, maf = c(pca_maf, 1), .progress = TRUE)

variant_filterlist <- read_tsv(variant_qc_summary_tsv)

message("Variant-level missing genotype rate filtering")
variant_filter_missingness <- variant_filterlist %$% variant.id
# variant.missing.rate <- missingGenotypeRate(gds, margin="by.variant")
# variant.filter <- seqGetData(gds, "variant.id")[variant.missing.rate <= snp_missing_rate]
seqSetFilter(gds, variant.id = variant_filter_missingness)

df.pca.metrics["Variant_missing_rate", ] <- c(snp_missing_rate, length(variant_filter_missingness), length(seqGetData(gds, "sample.id")))
# df.pca.metrics["Variant_missing_rate", ] <- c(snp_missing_rate, length(seqGetData(gds, "variant.id")), length(seqGetData(gds, "sample.id")))

message("Variant-level heterozygosity filtering")
variant_filter_heterozygosity <- variant_filterlist |>
                                      dplyr::filter(pass_heterozygosity == TRUE) %$% variant.id         
# variant.heterozygosity <- heterozygosity(gds, margin="by.variant")
# variant.filter <- seqGetData(gds, "variant.id")[variant.heterozygosity <= 0.6]
seqSetFilter(gds, variant.id = variant_filter_heterozygosity)

df.pca.metrics["Variant_heterozygosity_failed", ] <- c(">0.6", length(which(variant_filterlist$heterozygosity > 0.6)), 0)
# df.pca.metrics["Variant_heterozygosity_failed", ] <- c(">0.6", length(which(variant.heterozygosity > 0.6)), 0)

# message("Sample-level missing genotype rate filtering")
# sample.missing.rate <- missingGenotypeRate(gds, margin="by.sample")
# sample.filter <- seqGetData(gds, "sample.id")[sample.missing.rate <= snp_missing_rate]
# seqSetFilter(gds, sample.id = sample.filter)

# df.pca.metrics["Sample_missing_rate", ] <- c(snp_missing_rate, length(seqGetData(gds, "variant.id")), length(seqGetData(gds, "sample.id")))

df.pca.metrics |>
  rownames_to_column(var = "Parameter") |>
  write_tsv(pca_qc_summary_tsv)

df.pca.metrics

seqGDS2SNP(gds, pca_snp)

seqClose(gds)

sample_pheno %>% write_tsv(pca_smp)

sessionInfo()
