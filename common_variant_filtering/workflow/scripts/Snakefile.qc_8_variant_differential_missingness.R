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
# variant_missingness_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/test/adsp.r4.wgs_cc.all.chrall.variant_missingness_filter.tsv"
# variant_diff_missingness_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/adsp.r4.wgs_cc.all.variant_diff_missingness.tsv"
# qc_summary_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/qc_metrics.tsv"
# snp_maf <- 0.01
# missing.rate <- 0.1
# refp.super_pops <- "all"
# chr <- 21

gds_fn <- snakemake@input$gds
sample_pheno_tsv <- snakemake@input$sample_pheno_tsv
variant_missingness_tsv <- snakemake@input$variant_missingness_tsv
variant_diff_missingness_tsv <- snakemake@output$variant_diff_missingness_tsv
snp_maf <- snakemake@params$snp_maf
missing.rate <- snakemake@params$missing_rate
refp.super_pops <- snakemake@params$refp_super_pops
chr <- snakemake@wildcards$chr

## Script
`%notin%` <- Negate(`%in%`)

gds <- seqOpen(gds_fn)
# seqResetFilter(gds)

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

## Filter for samples in pheno file and biallelic variants
message("Filter for samples in pheno file and biallelic variants only")
seqSetFilter(gds, sample.id = sample_pheno$sample.id, variant.sel = isSNV(gds, biallelic=TRUE))

# seqSetFilterChrom(gds, include = chr, verbose = TRUE)

## Filter for variant missing rate
variant_metrics_filter <- read_tsv(variant_missingness_tsv) |>
  separate(variant.id, into = c("CHROM", "POS", "REF", "ALT"), sep = ":", remove = FALSE) |>
  dplyr::filter(CHROM == chr, !is.na(pass_missingness))

seqSetFilter(gds, variant.id = variant_metrics_filter$variant.id)

message("Calculate differential missingness by population and diagnosis")
variant.id <- seqGetData(gds, "variant.id")

## Get sample lists
get_sampleid_list <- function(sample_pheno, refp.super_pops) {
  sampleid_list <- list()
  
  # Always include all cases and controls
  sampleid_list[["cases_all"]] <- sample_pheno %>% 
    dplyr::filter(diagnosis == "AD") %$% sample.id
  sampleid_list[["controls_all"]] <- sample_pheno %>% 
    dplyr::filter(diagnosis == "Control") %$% sample.id
  
  # Define populations to filter
  pop_map <- list(
    "all" = c("EUR", "AFR", "AMR"),
    "EUR" = "EUR",
    "AFR" = "AFR",
    "AMR" = "AMR",
    "EUR_AMR" = c("EUR", "AMR")
  )
  
  pops_to_use <- pop_map[[refp.super_pops]]
  
  if (!is.null(pops_to_use)) {
    for (pop in pops_to_use) {
      sampleid_list[[paste0("cases_", pop)]] <- sample_pheno %>%
        dplyr::filter(admixture_super_pop_max == pop & diagnosis == "AD") %$% sample.id
      sampleid_list[[paste0("controls_", pop)]] <- sample_pheno %>%
        dplyr::filter(admixture_super_pop_max == pop & diagnosis == "Control") %$% sample.id
    }
  }
  
  return(sampleid_list)
}

sampleid_list <- get_sampleid_list(sample_pheno = sample_pheno, refp.super_pops = refp.super_pops)

## Get differential missingness
compute_missingness_by_pop <- function(gds, sampleid_list, refp.super_pops, variant_sel) {
  # Set variant filter
  seqSetFilter(gds, variant.id = variant_sel, verbose = TRUE)
  
  # Get current variant selection (after filtering)
  variant.id <- seqGetData(gds, "variant.id")
  
  # Always include these
  missingness <- list(
    variant.missing.rate = seqMissing(gds, per.variant = TRUE),
    variant.diff_missingness.rate.cases_all = applyMethod(
      gds, missingGenotypeRate, margin = "by.variant", use.names = TRUE,
      sample = sampleid_list$cases_all),
    variant.diff_missingness.rate.controls_all = applyMethod(
      gds, missingGenotypeRate, margin = "by.variant", use.names = TRUE,
      sample = sampleid_list$controls_all)
  )
  
  # Define populations per setting
  pop_map <- list(
    "all" = c("EUR", "AFR", "AMR"),
    "EUR" = "EUR",
    "AFR" = "AFR",
    "AMR" = "AMR",
    "EUR_AMR" = c("EUR", "AMR")
  )
  
  pops_to_include <- pop_map[[refp.super_pops]]
  
  # Add selected populations
  if (!is.null(pops_to_include)) {
    for (pop in pops_to_include) {
      missingness[[paste0("variant.diff_missingness.rate.cases_", pop)]] <- applyMethod(
        gds, missingGenotypeRate, margin = "by.variant", use.names = TRUE,
        variant = variant_sel,
        sample = sampleid_list[[paste0("variant.diff_missingness.rate.cases_", pop)]]
      )
      missingness[[paste0("variant.diff_missingness.rate.controls_", pop)]] <- applyMethod(
        gds, missingGenotypeRate, margin = "by.variant", use.names = TRUE,
        variant = variant_sel,
        sample = sampleid_list[[paste0("variant.diff_missingness.rate.controls_", pop)]]
      )
    }
  }
  
  # Combine into dataframe
  variant.diff_missingness.pop_dx <- cbind(variant.id = variant.id, as.data.frame(missingness))
  
  return(variant.diff_missingness.pop_dx)
}

variant.diff_missingness.pop_dx <- compute_missingness_by_pop(
  gds = gds,
  sampleid_list = sampleid_list,
  variant_sel = variant_metrics_filter$variant.id,
  refp.super_pops = refp.super_pops)

head(variant.diff_missingness.pop_dx)

## Compute variant missingness counts
compute_variant_missingness_counts <- function(variant_data, sampleid_list, refp.super_pops) {
  message("Variant-level missingness for all")
  result_list <- list()
  
  # Always compute "all"
  result_list[["all"]] <- variant_data %>%
    dplyr::select(variant.id, variant.missing.rate,
                  variant.diff_missingness.rate.cases_all,
                  variant.diff_missingness.rate.controls_all) %>%
    dplyr::mutate(
      missing_count.case = length(sampleid_list$cases_all) * as.numeric(variant.diff_missingness.rate.cases_all),
      present_count.case = length(sampleid_list$cases_all) - missing_count.case,
      missing_count.control = length(sampleid_list$controls_all) * as.numeric(variant.diff_missingness.rate.controls_all),
      present_count.control = length(sampleid_list$controls_all) - missing_count.control
    )
  
  # Define populations based on setting
  pop_map <- list(
    "all" = c("EUR", "AFR", "AMR"),
    "EUR" = "EUR",
    "AFR" = "AFR",
    "AMR" = "AMR",
    "EUR_AMR" = c("EUR", "AMR")
  )
  
  pops <- pop_map[[refp.super_pops]]
  if (!is.null(pops)) {
    for (pop in pops) {
      message(paste0("Variant-level missingness for ", pop))
      result_list[[pop]] <- variant_data %>%
        dplyr::select(variant.id, variant.missing.rate,
                      !!rlang::sym(paste0("variant.diff_missingness.rate.cases_", pop)),
                      !!rlang::sym(paste0("variant.diff_missingness.rate.controls_", pop))) %>%
        dplyr::rename(
          case_rate = !!rlang::sym(paste0("variant.diff_missingness.rate.cases_", pop)),
          control_rate = !!rlang::sym(paste0("variant.diff_missingness.rate.controls_", pop))
        ) %>%
        dplyr::mutate(
          missing_count.case = length(sampleid_list[[paste0("cases_", pop)]]) * as.numeric(case_rate),
          present_count.case = length(sampleid_list[[paste0("cases_", pop)]]) - missing_count.case,
          missing_count.control = length(sampleid_list[[paste0("controls_", pop)]]) * as.numeric(control_rate),
          present_count.control = length(sampleid_list[[paste0("controls_", pop)]]) - missing_count.control
        )
    }
  }
  
  return(result_list)
}

variant.diff_missingness_counts.pop_dx <- compute_variant_missingness_counts(
  variant_data = variant.diff_missingness.pop_dx, 
  sampleid_list = sampleid_list, 
  refp.super_pops = refp.super_pops)

message("Conduct Fisher Exact Test per variant")
## Function to run Fisher test per variant
fisher_test_variant_level <- function(variant.id, case_missing, control_missing, case_present, control_present) {
  # Ensure inputs are numeric and finite
  if (is.finite(case_missing) && case_missing >= 0 &&
      is.finite(control_missing) && control_missing >= 0 &&
      is.finite(case_present) && case_present >= 0 &&
      is.finite(control_present) && control_present >= 0) {
    
    # Round counts to the nearest integer to satisfy fisher.test requirements
    contingency_table <- matrix(
      round(c(case_missing, case_present,
              control_missing, control_present)),
      nrow = 2, byrow = TRUE,
      dimnames = list(c("case", "control"), c("missing", "present"))
    )
    
    test_results <- fisher.test(
      contingency_table,
      alternative = "two.sided",
      conf.int = TRUE,
      conf.level = 0.95
    )
    
    c(variant.id = variant.id, ci = test_results$conf.int, p = test_results$p.value,
      case_missing = contingency_table[1, 1],
      control_missing = contingency_table[2, 1],
      case_present = contingency_table[1, 2],
      control_present = contingency_table[2, 2])
  } else {
    c(variant.id = variant.id, ci = c(NA, NA), p = NA,
      case_missing = case_missing, control_missing = control_missing,
      case_present = case_present, control_present = control_present)
  }
}

## Fisher test on all variants in one population
fisher_test_allvars <- function(tbl, stratum) {
  tbl %>%
    dplyr::select(
      variant.id = variant.id,
      case_missing = missing_count.case,
      control_missing = missing_count.control,
      case_present = present_count.case,
      control_present = present_count.control
    ) %>%
    pmap_dfr(fisher_test_variant_level) %>%
    mutate(stratum = stratum)
}

## Define populations based on input
strata <- switch(
  refp.super_pops,
  "all" = c("EUR", "AFR", "AMR"),
  "EUR_AMR" = c("EUR", "AMR"),
  refp.super_pops  # default fallback
)

## Map the correct data frames from the variant.diff_missingness_counts.pop_dx list
fisher_input <- setNames(
  lapply(strata, function(pop) variant.diff_missingness_counts.pop_dx[[pop]]),
  strata
)

## Run Fisher test and pivot results
message("Get Fisher results for differential missingness")
fisher_results <- fisher_input %>%
  imap_dfr(fisher_test_allvars) %>%
  dplyr::select(variant.id, p, stratum) %>%
  pivot_wider(id_cols = variant.id, names_from = stratum, values_from = p)

head(fisher_results)

diff_miss_cols <- setNames(paste0("DIFF_MISS_", strata), strata)

message("Output variant differential missingness metrics file")
head(fisher_results)

variant_diff_missingness <- fisher_results %>%
  dplyr::rename_with(
    .fn = ~ paste0("DIFF_MISS_", .),
    .cols = all_of(strata[strata %in% names(fisher_results)])
  )

## Write output
head(variant_diff_missingness)

variant_diff_missingness %>% 
  readr::write_tsv(variant_diff_missingness_tsv)

seqClose(gds)

# showfile.gds(closeall=TRUE)

sessionInfo()
