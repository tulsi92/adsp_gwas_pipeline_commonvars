set.seed(19690329)

library(tidyverse)
library(magrittr)
library(GGally)
library(RColorBrewer)
library(SeqVarTools)
library(SNPRelate)
library(GENESIS)
library(readr)

# THREADS <- 40
# adsp_snp_gds <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/common_all/adsp.r4.wgs_cc.all.chrall_pca.snp"
# adsp_smp_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/common_all/adsp.r4.wgs_cc.all.chrall_pca.smp"
# refp_snp_gds <- "output/adsp.r4.wgs/cc.all/refp.gnomad.v3.1/all/refp.chrall_pca.snp"
# refp_smp_tsv <- "output/adsp.r4.wgs/cc.all/refp.gnomad.v3.1/all/refp.chrall_pca.smp"
# sample_pheno_tsv <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/common_all/adsp.r4.wgs_cc.all.sample_pheno.tsv"
# cov_mat_rds <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/common_all/adsp.r4.wgs_cc.all.cov_mat.rds"
# king_obj_rds <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/common_all/pca/adsp.r4.wgs_cc.all.king_obj.rds"
# pcair_obj_uncorrected_rds <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/common_all/pca/adsp.r4.wgs_cc.all.pcair_obj_uncorrected.rds"
# pcrelate_obj_uncorrected_rds <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/common_all/pca/adsp.r4.wgs_cc.all.pcrelate_obj_uncorrected.rds"
# pcair_obj_corrected_rds <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/common_all/pca/adsp.r4.wgs_cc.all.pcair_obj_corrected.rds"
# pcrelate_obj_corrected_rds <- "output/adsp.r4.wgs/cc.all/diagnosis_code/all/common_all/pca/adsp.r4.wgs_cc.all.pcrelate_obj_corrected.rds"
# pcs <- 10

THREADS <- snakemake@threads
adsp_snp_gds <- snakemake@input$adsp_snp
adsp_smp_tsv <- snakemake@input$adsp_smp
refp_snp_gds <- snakemake@input$refp_snp
refp_smp_tsv <- snakemake@input$refp_smp
sample_pheno_tsv <- snakemake@output$sample_pheno_tsv
cov_mat_rds <- snakemake@output$cov_mat_rds
king_obj_rds <- snakemake@output$king_obj_rds
pcair_obj_uncorrected_rds <- snakemake@output$pcair_obj_uncorrected_rds
pcrelate_obj_uncorrected_rds <- snakemake@output$pcrelate_obj_uncorrected_rds
pcair_obj_corrected_rds <- snakemake@output$pcair_obj_corrected_rds
pcrelate_obj_corrected_rds <- snakemake@output$pcrelate_obj_corrected_rds
pcs <- snakemake@params$pcs

## Edit PCs in config after manual inspection of plots
seqParallelSetup(cluster=THREADS)

## Read sample phenotype information for input and reference populations
adsp_smp <- read_tsv(adsp_smp_tsv) %>%
  mutate(population = case_when(population == "NHW" ~ "EUR",
                         population == "HI" ~ "AMR",
                         population == "AA" ~ "AFR",
                         TRUE ~ population),
         super_pop = "ADSP")

refp_smp <- read_tsv(refp_smp_tsv) %>%
  dplyr::rename(population = pop)

comb_tsv <- rbind(adsp_smp %>% dplyr::select(sample.id, super_pop, population),
                 refp_smp %>% dplyr::select(sample.id, super_pop, population))

message("Combining reference and sample genotypes")
comb_snp_gds <- tempfile()

snpgdsCombineGeno(c(adsp_snp_gds, refp_snp_gds), comb_snp_gds)

comb_snp <- snpgdsOpen(comb_snp_gds)

message("Pruning to get SNPs for PCAir/PCRelate")
## Sliding window 10e6 based on TOPMed; Error does not support multithreading anymore
comb_snp.pruned <- snpgdsLDpruning(comb_snp,
                                   method = "corr",
                                   slide.max.bp = 10e6,
                                   ld.threshold = sqrt(0.1))

comb_snp.pruned.snp.id <- unlist(comb_snp.pruned, use.names = FALSE)

### KING
## LD pruning not recommended in KING but TOPMed uses pruned set
## Do not prune or filter any "good" SNPs that pass QC prior to any KING inference, unless the number of variants
## is too large for computer memory, e.g. >100 million in WGS, in which case rare variants can be filtered out
message("Calculating kinship matrix using KING on unpruned data")
king.obj <- snpgdsIBDKING(comb_snp, num.thread = THREADS) #NEVER prune for KING

kinship <- snpgdsIBDSelection(king.obj)

saveRDS(king.obj, king_obj_rds)

king.mat <- kingToMatrix(king.obj)

## Convert SNPRelate GDS to SeqArray GDS
comb_gds_gds <- tempfile()

## Close GDS otherwise the next function call fails!
snpgdsClose(comb_snp)

seqSNP2GDS(comb_snp_gds, comb_gds_gds)

comb_gds <- seqOpen(comb_gds_gds)

### PCAir and PCRelate
message("Performing initial PCA")
pcair.obj <- pcair(comb_gds,
                  snp.include = comb_snp.pruned.snp.id, # Only use pruned SNPS. GOOD
                  kinobj = king.mat, #kin.thresh=2^(-9/2),
                  divobj = king.mat, #div.thresh=-2^(-9/2),
                  num.cores = THREADS)

summary(pcair.obj)

saveRDS(pcair.obj, pcair_obj_uncorrected_rds)

message("Extracting PCAir results")
sample_annot <- as_tibble(pcair.obj$vectors, rownames = NA) %>%
  rownames_to_column("sample.id") %>%
  left_join(comb_tsv, by = "sample.id") %>%
  mutate(super_pop = fct_relevel(super_pop, "ADSP"))

message("Running inital PCRelate")
comb_gds.iter <- SeqVarBlockIterator(comb_gds)

pcrelate.obj <- pcrelate(comb_gds.iter,
                        pcs = pcair.obj$vectors[, 1:pcs],
                        training.set = pcair.obj$unrels,
                        BPPARAM = BiocParallel::SerialParam())

saveRDS(pcrelate.obj, pcrelate_obj_uncorrected_rds)

pcrelate.mat <- pcrelateToMatrix(pcrelate.obj, scaleKin=1)

resetIterator(comb_gds.iter)

message("Rerunning PCAir with corrected PCRelate kinship matrix")
pcair.obj <- pcair(comb_gds,
                  snp.include = comb_snp.pruned.snp.id,
                  kinobj = pcrelate.mat, #kin.thresh=2^(-9/2),
                  divobj = king.mat, #div.thresh=-2^(-9/2),
                  num.cores = THREADS * 2)

saveRDS(pcair.obj, pcair_obj_corrected_rds)

summary(pcair.obj)

sample_annot <- as_tibble(pcair.obj$vectors, rownames = NA) %>% #, rownames = "sample.id", .name_repair = "check_unique") %>%
  rownames_to_column("sample.id") %>%
  left_join(comb_tsv, by = "sample.id") %>%
  mutate(super_pop = fct_relevel(super_pop, "ADSP"))

# stopifnot(sample_annot %>% dplyr::filter(!complete.cases(.)) %>% nrow == 0)

message("Rerunning PCRelate with corrected PCAir PCs")
pcrelate.obj <- pcrelate(comb_gds.iter,
                        pcs = pcair.obj$vectors[, 1:pcs],
                        training.set = pcair.obj$unrels)

message("STOP!")
message("Manually inspect PCA plots and edit number of analysis PCs in config appropriately")
adsp.sample.id <- sample_annot %>% filter(sample.id %in% adsp_smp$sample.id) %$% sample.id

saveRDS(pcrelate.obj, pcrelate_obj_corrected_rds)

## GRM
cov_mat <- pcrelateToMatrix(pcrelate.obj, scaleKin=2, sample.include = adsp.sample.id)

saveRDS(cov_mat, cov_mat_rds)

## Filter ADSP sample list to include only QC'ed samples
sample_annot <- adsp_smp %>% 
  dplyr::filter(sample.id %in% adsp.sample.id) %>% 
  left_join(as_tibble(pcair.obj$vectors, rownames = "sample.id"), by = "sample.id")

sample_annot

sample_annot %>% 
  write_tsv(sample_pheno_tsv)

seqClose(comb_gds)

unlink(comb_gds)

showfile.gds(closeall=TRUE)