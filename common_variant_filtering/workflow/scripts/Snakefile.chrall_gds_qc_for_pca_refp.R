set.seed(19620321)

library(tidyverse)
library(magrittr)
library(SeqVarTools)
library(SNPRelate)
library(readr)

# gds_fn <- "/sc/arion/projects/load/users/patelt16/projects/adsp/apoe_extremes/070822.ref.vcf_to_gds/output/refp.gnomad.v3.1.chrall.gds"
# sample_pheno_tsv <- "/sc/arion/projects/load/users/patelt16/projects/adsp/apoe_extremes/070822.ref.vcf_to_gds/input/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta_populations.tsv"
# variant_qc_summary_tsv <- ""
# pca_snp <- "output/adsp_pca.snp"
# pca_smp <- "output/adsp_pca.smp"
# pca_maf <- 0.05
# snp_missing_rate <- 0.05
# refp.super_pops <- "all"

gds_fn <- snakemake@input$gds
sample_pheno_tsv <- snakemake@input$sample_pheno_tsv
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

seqSetFilterCond(gds, maf = c(pca_maf, 1), .progress = TRUE)

message("Variant-level missing genotype rate filtering")
variant.missing.rate <- missingGenotypeRate(gds, margin="by.variant")
variant.filter <- seqGetData(gds, "variant.id")[variant.missing.rate <= snp_missing_rate]
seqSetFilter(gds, variant.id = variant.filter)

df.pca.metrics["Variant_missing_rate", ] <- c(snp_missing_rate, length(seqGetData(gds, "variant.id")), length(seqGetData(gds, "sample.id")))

message("Variant-level heterozygosity filtering")
variant.heterozygosity <- heterozygosity(gds, margin="by.variant")
variant.filter <- seqGetData(gds, "variant.id")[variant.heterozygosity <= 0.6]
seqSetFilter(gds, variant.id = variant.filter)

df.pca.metrics["Variant_heterozygosity_failed", ] <- c(">0.6", length(which(variant.heterozygosity > 0.6)), 0)

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
