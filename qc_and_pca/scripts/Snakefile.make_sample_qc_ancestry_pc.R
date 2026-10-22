
#!/usr/bin/Rscript

##load packages and set variables

library(readr)

library(dplyr)

library(stringr)

args <- commandArgs(trailingOnly = TRUE)

adsp_sample_qc_name <- args[1]

adsp_ancestry_name <- args[2]

adsp_pc_name <- args[3]

adsp_sample_qc_mis_het_titv_inbreed_name <- args[4]

adsp_ancestry_clean_name <- args[5]

adsp_pc_clean_name <- args[6]

imiss_max <- as.numeric(args[7])



##import and process adsp sample qc data

adsp_sample_qc <- read_tsv(adsp_sample_qc_name, guess_max = 1e7) %>%

    rename(IID = sample.id) %>%

    mutate(FID = IID) %>%

    mutate(pass_missing = ifelse(missing_rate < imiss_max, TRUE, FALSE))

adsp_sample_qc_mis <- adsp_sample_qc %>%

    filter(pass_missing == FALSE) %>%

    mutate(why = "mis") %>%

    select(FID, IID, why)

adsp_sample_qc_het <- adsp_sample_qc %>%

    filter(pass_hethom == FALSE) %>%

    mutate(why = "het") %>%

    select(FID, IID, why)

adsp_sample_qc_titv <- adsp_sample_qc %>%

    filter(pass_titv == FALSE) %>%

    mutate(why = "titv") %>%

    select(FID, IID, why)

adsp_sample_qc_inbreed <- adsp_sample_qc %>%

    filter(pass_inbreed_coeff == FALSE) %>%

    mutate(why = "inbreed") %>%

    select(FID, IID, why)

adsp_sample_qc_mis_het_titv_inbreed <- bind_rows(adsp_sample_qc_mis, adsp_sample_qc_het, adsp_sample_qc_titv, adsp_sample_qc_inbreed) %>%

    arrange(FID, IID)



##import and process adsp ancestry data

adsp_ancestry <- read_tsv(adsp_ancestry_name, guess_max = 1e7)

adsp_ancestry_clean <- adsp_ancestry %>%

    mutate(FID = IID) %>%

    relocate(FID, IID) %>%

    arrange(FID, IID)



##import and process adsp pc data

adsp_pc <- read_tsv(adsp_pc_name, guess_max = 1e7)

adsp_pc_clean <- adsp_pc %>%

    rename(IID = sample.id) %>%

    mutate(FID = IID) %>%

    select(FID, IID, c(V1:V32)) %>%

    arrange(FID, IID)



##output adsp data

write_tsv(adsp_sample_qc_mis_het_titv_inbreed, adsp_sample_qc_mis_het_titv_inbreed_name)

write_tsv(adsp_ancestry_clean, adsp_ancestry_clean_name)

write_tsv(adsp_pc_clean, adsp_pc_clean_name, col_names = FALSE)

