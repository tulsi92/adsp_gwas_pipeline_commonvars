set.seed(19690329)

library(tidyverse)
library(magrittr)
library(SeqVarTools)
library(SNPRelate)
library(readr)

# gds_fn <- "output/adsp.r3.wgs/apoe33.extremes.nhw/diagnosis_code/EUR/adsp.r3.wgs_apoe33.extremes.nhw.chrall.filtered.gds"
# tsv <- "input/20230906.ADSP_Umbrella_ng00067.v8.wgs.nhw.apoe33.extremes.1063.tsv"
# gds_filtered_qc <- "output/adsp.gds"
# smp_filtered_qc <- "output/adsp.smp"
# snp_maf <- 0.01
# missing.rate <- 0.1
# refp.super_pops <- "EUR"
# qc_metrics_tsv <- "output/qc_metrics.tsv"
# qc_report_pdf <- "output/qc_report.pdf"

gds_fn <- snakemake@input$gds
tsv <- snakemake@input$tsv
gds_filtered_qc <- snakemake@output$gds
smp_filtered_qc <- snakemake@output$smp
variant_qc_flags_tsv <- snakemake@output$variant_qc_flags_tsv
qc_summary_tsv <- snakemake@output$qc_summary_tsv
qc_report_pdf <- snakemake@output$qc_report_pdf
snp_maf <- snakemake@params$snp_maf
missing.rate <- snakemake@params$missing_rate
refp.super_pops <- snakemake@params$refp_super_pops

`%notin%` <- Negate(`%in%`)

gds <- seqOpen(gds_fn)

tsv <- read_tsv(tsv)

if (length(refp.super_pops) > 0 & 'super_pop' %notin% colnames(tsv)) {
  tsv <- tsv %>% mutate(super_pop = case_when(population == "NHW" ~ "EUR",
                                             population == "AA" ~ "AFR",
                                             population == "HI" ~ "AMR",
                                             is.na(population) ~ NA,
                                             TRUE ~ "OTHER"))
}

if (length(refp.super_pops) > 0 & refp.super_pops != "all" & 'super_pop' %in% colnames(tsv)) {
  refp.super_pops_all <- unlist(strsplit(refp.super_pops, split = "_"))
  tsv %<>% filter(super_pop %in% c(refp.super_pops_all))
}

seqResetFilter(gds)

df.metrics <- data.frame(Threshold = numeric(0), N_Variants = numeric(0), N_Samples = numeric(0))
df.metrics["Total", ] <- c(Threshold = "null", N_Variants = length(seqGetData(gds, "variant.id")), N_Samples = length(seqGetData(gds, "sample.id")))

message("Filter for samples in pheno file and biallelic variants only")
seqSetFilter(gds, sample.id = tsv$sample.id, variant.sel = isSNV(gds, biallelic=TRUE))

df.metrics["Biallelic_and_analysis_samples", ] <- c("null", length(seqGetData(gds, "variant.id")), length(seqGetData(gds, "sample.id")))

pdf(qc_report_pdf, width = 1280, height = 960)
message("Variant-level missing genotype rate filtering")
variant.missing.rate <- missingGenotypeRate(gds, margin="by.variant")
hist(variant.missing.rate, breaks = 100)
abline(v=missing.rate, col = "red")
variant.filter_missing <- seqGetData(gds, "variant.id")[variant.missing.rate <= missing.rate]
seqSetFilter(gds, variant.id = variant.filter_missing)

df.metrics["Variant_missing_rate", ] <- c(missing.rate, length(seqGetData(gds, "variant.id")), length(seqGetData(gds, "sample.id")))

message("Sample-level missing genotype rate filtering")
sample.missing.rate <- missingGenotypeRate(gds, margin="by.sample")
hist(sample.missing.rate, breaks = 100)
abline(v=missing.rate, col = "red")
sample.filter_missing <- seqGetData(gds, "sample.id")[sample.missing.rate <= missing.rate]
seqSetFilter(gds, sample.id = sample.filter_missing)

df.metrics["Sample_missing_rate", ] <- c(missing.rate, length(seqGetData(gds, "variant.id")), length(seqGetData(gds, "sample.id")))

message("Sample-level hetero/homozygosity filtering")
sample.hethom <- hethom(gds)
hist(sample.hethom, breaks = 100)
sample.filter_hethom <- seqGetData(gds, "sample.id")[abs(scale(sample.hethom)[,1]) <= 3]
seqSetFilter(gds, sample.id = sample.filter_hethom)

df.metrics["Sample_hethom_rate", ] <- c("<=3", length(seqGetData(gds, "variant.id")), length(seqGetData(gds, "sample.id")))

message("Ti/Tv ratio should be approximately 2.1 for whole genome sequencing and 2.8 for whole exome sequencing!")
sample.titv <- titv(gds, by.sample=TRUE)
hist(sample.titv)
abline(v=c(2.1, 2.8), col = "red")
sample.filter_titv <- seqGetData(gds, "sample.id")[abs(scale(sample.titv)[,1]) <= 3]
seqSetFilter(gds, sample.id = sample.filter_titv)

df.metrics["Sample_titv_ratio", ] <- c("<=3", length(seqGetData(gds, "variant.id")), length(seqGetData(gds, "sample.id")))

message("Sample-level inbreeding coefficient filtering")
sample.inbreedCoeff <- inbreedCoeff(gds, margin="by.sample")
hist(sample.inbreedCoeff)
sample.filter_inbreeding <- seqGetData(gds, "sample.id")[abs(scale(sample.inbreedCoeff)[,1]) <= 3]
seqSetFilter(gds, sample.id = sample.filter_inbreeding)

df.metrics["Sample_inbreeding_coeff", ] <- c("<=3", length(seqGetData(gds, "variant.id")), length(seqGetData(gds, "sample.id")))

## calculate metrics but do not filter
message("Variant-level heterozygosity filtering")
variant.heterozygosity <- heterozygosity(gds, margin="by.variant")
hist(variant.heterozygosity, breaks = 100)
abline(v=0.6, col = "red")
variant.filter_heterozygosity <- seqGetData(gds, "variant.id")[variant.heterozygosity <= 0.6]
# seqSetFilter(gds, variant.id = variant.filter_heterozygosity)

head(df.metrics)

df.variant_metrics <- data.frame()
df.variant_metrics[seqGetData(gds, "variant.id"), "variant_heterozygosity"] <- variant.heterozygosity

df.metrics["Variant_heterozygosity_failed", ] <- c(">0.6", length(which(variant.heterozygosity > 0.6)), 0)

message("Variant-level Hardy-Weinberg Equilibrium (HWE)/inbreeding coefficient filtering")
variant.hwe <- hwe(gds)
observed <- -log10(sort(variant.hwe$p))
variant.hwe.perm <- hwe(gds, permute=TRUE)
expected <- -log10(sort(variant.hwe.perm$p))
plot(x=expected, y=observed, xlab="-log10(expected P)", ylab="-log10(observed P)")
abline(0, 1, col="red")
## hwe$f <- inbreeding coefficient <- 1 - observed heterozygosity / expected heterozygosity
hist(variant.hwe$f, breaks = 100)
variant.filter_hwe <- seqGetData(gds, "variant.id")[abs(scale(variant.hwe$f)[,1]) <= 3]
# seqSetFilter(gds, variant.id = variant.filter_hwe, action = "intersect")

df.variant_metrics[, "variant_hwe_overall"] <- variant.hwe$p
df.variant_metrics[, "MAF"] <- variant.hwe$afreq

df.metrics["Variant_HWE_failed", ] <- c("<0.05", length(which(variant.hwe$p < 0.05)), 0)

df.metrics %>% 
  rownames_to_column(var = "Parameter") %>% 
  write_tsv(qc_summary_tsv)

# grid.table(df.metrics)

dev.off()

tsv %<>% 
  filter(sample.id %in% seqGetData(gds, "sample.id")) %>% 
  write_tsv(smp_filtered_qc)

message("Differential missingness by population and diagnosis")
variant.id.ann <- seqGetData(gds, "variant.id")

sampleid_list <- list()

for (i in nrow(tsv)){
  print(i)
  if(refp.super_pops == "all") {
    sampleid_list[["cases_all"]] <- tsv %>% filter(diagnosis == "AD") %$% sample.id
    sampleid_list[["controls_all"]] <- tsv %>% filter(diagnosis == "Control") %$% sample.id
    sampleid_list[["cases_NHW"]] <- tsv %>% filter(population == "NHW" & diagnosis == "AD") %$% sample.id
    sampleid_list[["controls_NHW"]] <- tsv %>% filter(population == "NHW" & diagnosis == "Control") %$% sample.id
    sampleid_list[["cases_AA"]] <- tsv %>% filter(population == "AA" & diagnosis == "AD") %$% sample.id
    sampleid_list[["controls_AA"]] <- tsv %>% filter(population == "AA" & diagnosis == "Control") %$% sample.id
    sampleid_list[["cases_HI"]] <- tsv %>% filter(population == "HI" & diagnosis == "AD") %$% sample.id
    sampleid_list[["controls_HI"]] <- tsv %>% filter(population == "HI" & diagnosis == "Control") %$% sample.id
    }
  else if(refp.super_pops == "EUR") {
    sampleid_list[["cases_all"]] <- tsv %>% filter(diagnosis == "AD") %$% sample.id
    sampleid_list[["controls_all"]] <- tsv %>% filter(diagnosis == "Control") %$% sample.id
    sampleid_list[["cases_NHW"]] <- tsv %>% filter(population == "NHW" & diagnosis == "AD") %$% sample.id
    sampleid_list[["controls_NHW"]] <- tsv %>% filter(population == "NHW" & diagnosis == "Control") %$% sample.id
  } else if(refp.super_pops == "AFR") {
    sampleid_list[["cases_all"]] <- tsv %>% filter(diagnosis == "AD") %$% sample.id
    sampleid_list[["controls_all"]] <- tsv %>% filter(diagnosis == "Control") %$% sample.id
    sampleid_list[["cases_AA"]] <- tsv %>% filter(population == "AA" & diagnosis == "AD") %$% sample.id
    sampleid_list[["controls_AA"]] <- tsv %>% filter(population == "AA" & diagnosis == "Control") %$% sample.id
  } else if(refp.super_pops == "AMR") {
    sampleid_list[["cases_all"]] <- tsv %>% filter(diagnosis == "AD") %$% sample.id
    sampleid_list[["controls_all"]] <- tsv %>% filter(diagnosis == "Control") %$% sample.id
    sampleid_list[["cases_HI"]] <- tsv %>% filter(population == "HI" & diagnosis == "AD") %$% sample.id
    sampleid_list[["controls_HI"]] <- tsv %>% filter(population == "HI" & diagnosis == "Control") %$% sample.id
  } else if(refp.super_pops == "EUR_AMR") {
    sampleid_list[["cases_all"]] <- tsv %>% filter(diagnosis == "AD") %$% sample.id
    sampleid_list[["controls_all"]] <- tsv %>% filter(diagnosis == "Control") %$% sample.id
    sampleid_list[["cases_NHW"]] <- tsv %>% filter(population == "NHW" & diagnosis == "AD") %$% sample.id
    sampleid_list[["controls_NHW"]] <- tsv %>% filter(population == "NHW" & diagnosis == "Control") %$% sample.id
    sampleid_list[["cases_HI"]] <- tsv %>% filter(population == "HI" & diagnosis == "AD") %$% sample.id
    sampleid_list[["controls_HI"]] <- tsv %>% filter(population == "HI" & diagnosis == "Control") %$% sample.id
  }
}

if(refp.super_pops == "all") {
  variant.diff_missingness.rate.cases_all <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$cases_all)
  variant.diff_missingness.rate.controls_all <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$controls_all)
  variant.diff_missingness.rate.cases_NHW <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$cases_NHW)
  variant.diff_missingness.rate.controls_NHW <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$controls_NHW)
  variant.diff_missingness.rate.cases_AA <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$cases_AA)
  variant.diff_missingness.rate.controls_AA <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$controls_AA)
  variant.diff_missingness.rate.cases_HI <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$cases_HI)
  variant.diff_missingness.rate.controls_HI <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$controls_HI)
} else if(refp.super_pops == "EUR") {
    variant.diff_missingness.rate.cases_all <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$cases_all)
    variant.diff_missingness.rate.controls_all <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$controls_all)
    variant.diff_missingness.rate.cases_NHW <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$cases_NHW)
    variant.diff_missingness.rate.controls_NHW <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$controls_NHW)
} else if(refp.super_pops == "AFR") {
    variant.diff_missingness.rate.cases_all <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$cases_all)
    variant.diff_missingness.rate.controls_all <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$controls_all)
    variant.diff_missingness.rate.cases_AA <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$cases_AA)
    variant.diff_missingness.rate.controls_AA <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$controls_AA)
} else if(refp.super_pops == "AMR") {
    variant.diff_missingness.rate.cases_all <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$cases_all)
    variant.diff_missingness.rate.controls_all <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$controls_all)
    variant.diff_missingness.rate.cases_HI <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$cases_HI)
    variant.diff_missingness.rate.controls_HI <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$controls_HI)
} else if(refp.super_pops == "EUR_AMR") {
    variant.diff_missingness.rate.cases_all <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$cases_all)
    variant.diff_missingness.rate.controls_all <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$controls_all)
    variant.diff_missingness.rate.cases_NHW <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$cases_NHW)
    variant.diff_missingness.rate.controls_NHW <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$controls_NHW)
    variant.diff_missingness.rate.cases_HI <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$cases_HI)
    variant.diff_missingness.rate.controls_HI <- applyMethod(gds, missingGenotypeRate, margin="by.variant", use.names=TRUE, sample=sampleid_list$controls_HI)
}

if(refp.super_pops == "all") {
  variant.diff_missingness.pop_dx <- cbind(variant.id.ann, variant.missing.rate = seqMissing(gds, per.variant = TRUE), 
                                         variant.diff_missingness.rate.cases_all, variant.diff_missingness.rate.controls_all,
                                         variant.diff_missingness.rate.cases_NHW, variant.diff_missingness.rate.controls_NHW,
                                         variant.diff_missingness.rate.cases_AA, variant.diff_missingness.rate.controls_AA,
                                         variant.diff_missingness.rate.cases_HI, variant.diff_missingness.rate.controls_HI) %>% as.data.frame()
} else if(refp.super_pops == "EUR") {
  variant.diff_missingness.pop_dx <- cbind(variant.id.ann, variant.missing.rate = seqMissing(gds, per.variant = TRUE), 
                                         variant.diff_missingness.rate.cases_all, variant.diff_missingness.rate.controls_all,
                                         variant.diff_missingness.rate.cases_NHW, variant.diff_missingness.rate.controls_NHW) %>% as.data.frame()
} else if(refp.super_pops == "AFR") {
  variant.diff_missingness.pop_dx <- cbind(variant.id.ann, variant.missing.rate = seqMissing(gds, per.variant = TRUE), 
                                         variant.diff_missingness.rate.cases_all, variant.diff_missingness.rate.controls_all,
                                         variant.diff_missingness.rate.cases_AA, variant.diff_missingness.rate.controls_AA) %>% as.data.frame()
} else if(refp.super_pops == "AMR") {
  variant.diff_missingness.pop_dx <- cbind(variant.id.ann, variant.missing.rate = seqMissing(gds, per.variant = TRUE), 
                                         variant.diff_missingness.rate.cases_all, variant.diff_missingness.rate.controls_all,
                                         variant.diff_missingness.rate.cases_HI, variant.diff_missingness.rate.controls_HI) %>% as.data.frame()
} else if(refp.super_pops == "EUR_AMR") {
  variant.diff_missingness.pop_dx <- cbind(variant.id.ann, variant.missing.rate = seqMissing(gds, per.variant = TRUE), 
                                         variant.diff_missingness.rate.cases_all, variant.diff_missingness.rate.controls_all,
                                         variant.diff_missingness.rate.cases_NHW, variant.diff_missingness.rate.controls_NHW,
                                         variant.diff_missingness.rate.cases_HI, variant.diff_missingness.rate.controls_HI) %>% as.data.frame()
}

message("Variant-level missingness for all")
variant.diff_missingness.rate.all.DF <- variant.diff_missingness.pop_dx %>%
  dplyr::select(variant.id.ann, variant.missing.rate, variant.diff_missingness.rate.cases_all, variant.diff_missingness.rate.controls_all) %>%
  mutate(missing_count.case = (length(sampleid_list$cases_all) * as.numeric(variant.diff_missingness.rate.cases_all)),
         present_count.case = (length(sampleid_list$cases_all) - missing_count.case),
         missing_count.control = length(sampleid_list$controls_all) * as.numeric(variant.diff_missingness.rate.controls_all),
         present_count.control = length(sampleid_list$controls_all) - missing_count.control)

if(refp.super_pops == "all") {
  message("Variant-level missingness for NHW")
  variant.diff_missingness.rate.NHW.DF <- variant.diff_missingness.pop_dx %>% 
    dplyr::select(variant.id.ann, variant.missing.rate, variant.diff_missingness.rate.cases_NHW, variant.diff_missingness.rate.controls_NHW) %>%
    mutate(missing_count.case = (length(sampleid_list$cases_NHW) * as.numeric(variant.diff_missingness.rate.cases_NHW)),
          present_count.case = (length(sampleid_list$cases_NHW) - missing_count.case), 
          missing_count.control = length(sampleid_list$controls_NHW) * as.numeric(variant.diff_missingness.rate.controls_NHW),
          present_count.control = length(sampleid_list$controls_NHW) - missing_count.control)
  message("Variant-level missingness for AA")
  variant.diff_missingness.rate.AA.DF <- variant.diff_missingness.pop_dx %>% 
    dplyr::select(variant.id.ann, variant.missing.rate, variant.diff_missingness.rate.cases_AA, variant.diff_missingness.rate.controls_AA) %>%
    mutate(missing_count.case = (length(sampleid_list$cases_AA) * as.numeric(variant.diff_missingness.rate.cases_AA)),
          present_count.case = (length(sampleid_list$cases_AA) - missing_count.case), 
          missing_count.control = length(sampleid_list$controls_AA) * as.numeric(variant.diff_missingness.rate.controls_AA),
          present_count.control = length(sampleid_list$controls_AA) - missing_count.control)
  message("Variant-level missingness for HI")
  variant.diff_missingness.rate.HI.DF <- variant.diff_missingness.pop_dx %>% 
    dplyr::select(variant.id.ann, variant.missing.rate, variant.diff_missingness.rate.cases_HI, variant.diff_missingness.rate.controls_HI) %>%
    mutate(missing_count.case = (length(sampleid_list$cases_HI) * as.numeric(variant.diff_missingness.rate.cases_HI)),
          present_count.case = (length(sampleid_list$cases_HI) - missing_count.case), 
          missing_count.control = length(sampleid_list$controls_HI) * as.numeric(variant.diff_missingness.rate.controls_HI),
          present_count.control = length(sampleid_list$controls_HI) - missing_count.control)
} else if(refp.super_pops == "EUR") {
  message("Variant-level missingness for NHW")
  variant.diff_missingness.rate.NHW.DF <- variant.diff_missingness.pop_dx %>% 
    dplyr::select(variant.id.ann, variant.missing.rate, variant.diff_missingness.rate.cases_NHW, variant.diff_missingness.rate.controls_NHW) %>%
    mutate(missing_count.case = (length(sampleid_list$cases_NHW) * as.numeric(variant.diff_missingness.rate.cases_NHW)),
          present_count.case = (length(sampleid_list$cases_NHW) - missing_count.case), 
          missing_count.control = length(sampleid_list$controls_NHW) * as.numeric(variant.diff_missingness.rate.controls_NHW),
          present_count.control = length(sampleid_list$controls_NHW) - missing_count.control)
} else if(refp.super_pops == "AFR") {
  message("Variant-level missingness for AA")
  variant.diff_missingness.rate.AA.DF <- variant.diff_missingness.pop_dx %>% 
    dplyr::select(variant.id.ann, variant.missing.rate, variant.diff_missingness.rate.cases_AA, variant.diff_missingness.rate.controls_AA) %>%
    mutate(missing_count.case = (length(sampleid_list$cases_AA) * as.numeric(variant.diff_missingness.rate.cases_AA)),
          present_count.case = (length(sampleid_list$cases_AA) - missing_count.case), 
          missing_count.control = length(sampleid_list$controls_AA) * as.numeric(variant.diff_missingness.rate.controls_AA),
          present_count.control = length(sampleid_list$controls_AA) - missing_count.control)
} else if(refp.super_pops == "AMR") {
  message("Variant-level missingness for HI")
  variant.diff_missingness.rate.HI.DF <- variant.diff_missingness.pop_dx %>% 
    dplyr::select(variant.id.ann, variant.missing.rate, variant.diff_missingness.rate.cases_HI, variant.diff_missingness.rate.controls_HI) %>%
    mutate(missing_count.case = (length(sampleid_list$cases_HI) * as.numeric(variant.diff_missingness.rate.cases_HI)),
          present_count.case = (length(sampleid_list$cases_HI) - missing_count.case), 
          missing_count.control = length(sampleid_list$controls_HI) * as.numeric(variant.diff_missingness.rate.controls_HI),
          present_count.control = length(sampleid_list$controls_HI) - missing_count.control)
} else if(refp.super_pops == "EUR_AMR") {
  message("Variant-level missingness for NHW")
  variant.diff_missingness.rate.NHW.DF <- variant.diff_missingness.pop_dx %>% 
    dplyr::select(variant.id.ann, variant.missing.rate, variant.diff_missingness.rate.cases_NHW, variant.diff_missingness.rate.controls_NHW) %>%
    mutate(missing_count.case = (length(sampleid_list$cases_NHW) * as.numeric(variant.diff_missingness.rate.cases_NHW)),
          present_count.case = (length(sampleid_list$cases_NHW) - missing_count.case), 
          missing_count.control = length(sampleid_list$controls_NHW) * as.numeric(variant.diff_missingness.rate.controls_NHW),
          present_count.control = length(sampleid_list$controls_NHW) - missing_count.control)
  message("Variant-level missingness for HI")
  variant.diff_missingness.rate.HI.DF <- variant.diff_missingness.pop_dx %>% 
    dplyr::select(variant.id.ann, variant.missing.rate, variant.diff_missingness.rate.cases_HI, variant.diff_missingness.rate.controls_HI) %>%
    mutate(missing_count.case = (length(sampleid_list$cases_HI) * as.numeric(variant.diff_missingness.rate.cases_HI)),
          present_count.case = (length(sampleid_list$cases_HI) - missing_count.case), 
          missing_count.control = length(sampleid_list$controls_HI) * as.numeric(variant.diff_missingness.rate.controls_HI),
          present_count.control = length(sampleid_list$controls_HI) - missing_count.control)
}

message("Conduct Fisher Exact Test per variant")
# fisher_test_variant_level <- function(var.id, case_missing, control_missing, case_present, control_present) {
#   contingency_table <- matrix(
#     c(case_missing, case_present,
#       control_missing, control_present),
#     dimnames = list(c("case", "control"),
#                     c("missing", "present")),
#     nrow = 2, byrow = TRUE)
#   test_results <- fisher.test(
#     contingency_table,
#     alternative = "two.sided",
#     conf.int = TRUE,
#     conf.level = 0.95)
#   c(var.id = var.id, ci = test_results$conf.int, p = test_results$p.value,
#     case_missing = case_missing, control_missing = control_missing,
#     case_present = case_present, control_present = control_present)
# }

## Updated code to catch infinite error
fisher_test_variant_level <- function(var.id, case_missing, control_missing, case_present, control_present) {
  if (is.finite(case_missing) && case_missing >= 0 &&
      is.finite(control_missing) && control_missing >= 0 &&
      is.finite(case_present) && case_present >= 0 &&
      is.finite(control_present) && control_present >= 0) {
    contingency_table <- matrix(
      c(case_missing, case_present,
        control_missing, control_present),
      dimnames = list(c("case", "control"),
                      c("missing", "present")),
      nrow = 2, byrow = TRUE)
    test_results <- fisher.test(
      contingency_table,
      alternative = "two.sided",
      conf.int = TRUE,
      conf.level = 0.95)
    c(var.id = var.id, ci = test_results$conf.int, p = test_results$p.value,
      case_missing = case_missing, control_missing = control_missing,
      case_present = case_present, control_present = control_present)
  } else {
    c(var.id = var.id, ci = c(NA, NA), p = NA,
      case_missing = case_missing, control_missing = control_missing,
      case_present = case_present, control_present = control_present)
  }
}

fisher_test_allvars <- function(tbl, stratum) {
  tbl %>%
    dplyr::select(var.id = variant.id.ann,
                  case_missing = missing_count.case,
                  control_missing = missing_count.control,
                  case_present = present_count.case,
                  control_present = present_count.control) %>%
    pmap_dfr(fisher_test_variant_level) %>%
    mutate(stratum = stratum)
}

if(refp.super_pops == "all") {
  message("Get Fisher results for differential missingness")
  fisher_results <- list(NHW = variant.diff_missingness.rate.NHW.DF,
                        AA = variant.diff_missingness.rate.AA.DF,
                        HI = variant.diff_missingness.rate.HI.DF) %>%
    imap_dfr(fisher_test_allvars) %>%
    dplyr::select(var.id, p, stratum) %>%
    pivot_wider(id_cols = var.id, names_from = stratum, values_from = p)

  head(fisher_results)

  message("Output variant QC metrics file")
  head(df.variant_metrics)

  df.variant_metrics_all <- fisher_results %>% dplyr::select(var.id, DIFF_MISS_NHW = NHW, DIFF_MISS_AA = AA, DIFF_MISS_HI = HI) %>%
    left_join(df.variant_metrics %>% rownames_to_column(var = "var.id"), by = "var.id")
} else if(refp.super_pops == "EUR") {
  message("Get Fisher results for differential missingness")
  fisher_results <- list(NHW = variant.diff_missingness.rate.NHW.DF) %>%
    imap_dfr(fisher_test_allvars) %>%
    dplyr::select(var.id, p, stratum) %>%
    pivot_wider(id_cols = var.id, names_from = stratum, values_from = p)

  head(fisher_results)

  message("Output variant QC metrics file")
  head(df.variant_metrics)

  df.variant_metrics_all <- fisher_results %>% dplyr::select(var.id, DIFF_MISS_NHW = NHW) %>%
    left_join(df.variant_metrics %>% rownames_to_column(var = "var.id"), by = "var.id")
} else if(refp.super_pops == "AFR") {
  message("Get Fisher results for differential missingness")
  fisher_results <- list(AA = variant.diff_missingness.rate.AA.DF) %>%
    imap_dfr(fisher_test_allvars) %>%
    dplyr::select(var.id, p, stratum) %>%
    pivot_wider(id_cols = var.id, names_from = stratum, values_from = p)

  head(fisher_results)

  message("Output variant QC metrics file")
  head(df.variant_metrics)

  df.variant_metrics_all <- fisher_results %>% dplyr::select(var.id, DIFF_MISS_AA = AA) %>%
    left_join(df.variant_metrics %>% rownames_to_column(var = "var.id"), by = "var.id")
} else if(refp.super_pops == "AMR") {
  message("Get Fisher results for differential missingness")
  fisher_results <- list(HI = variant.diff_missingness.rate.HI.DF) %>%
    imap_dfr(fisher_test_allvars) %>%
    dplyr::select(var.id, p, stratum) %>%
    pivot_wider(id_cols = var.id, names_from = stratum, values_from = p)

  head(fisher_results)

  message("Output variant QC metrics file")
  head(df.variant_metrics)

  df.variant_metrics_all <- fisher_results %>% dplyr::select(var.id, DIFF_MISS_HI = HI) %>%
    left_join(df.variant_metrics %>% rownames_to_column(var = "var.id"), by = "var.id")
} else if(refp.super_pops == "EUR_AMR") {
    message("Get Fisher results for differential missingness")
  fisher_results <- list(NHW = variant.diff_missingness.rate.NHW.DF,
                        HI = variant.diff_missingness.rate.HI.DF) %>%
    imap_dfr(fisher_test_allvars) %>%
    dplyr::select(var.id, p, stratum) %>%
    pivot_wider(id_cols = var.id, names_from = stratum, values_from = p)

  head(fisher_results)

  message("Output variant QC metrics file")
  head(df.variant_metrics)

  df.variant_metrics_all <- fisher_results %>% dplyr::select(var.id, DIFF_MISS_NHW = NHW, DIFF_MISS_HI = HI) %>%
    left_join(df.variant_metrics %>% rownames_to_column(var = "var.id"), by = "var.id")
}

head(df.variant_metrics_all)

df.variant_metrics_all %>%
  write_tsv(variant_qc_flags_tsv)

message("Save QC filtered GDS")
seqExport(gds, gds_filtered_qc, verbose = TRUE)

seqClose(gds)

# showfile.gds(closeall=TRUE)

sessionInfo()
