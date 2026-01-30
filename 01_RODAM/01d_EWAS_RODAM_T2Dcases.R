############################################
# RODAM EWAS Adiponectin
# 01d_EWAS_RODAM_T2Dcases.R
# Author: Muhulo Muhau Mungamba
# Date: 2025-03-16
#
# EWAS of log-adiponectin in RODAM
# T2D cases only (n = 112)
# Fully adjusted model
############################################

rm(list = ls())
gc()

############################
# Load required packages
############################

library(IlluminaHumanMethylation450kmanifest)
library(minfi)
library(minfiData)
library(limma)
library(affy)
library(dplyr)

options(stringsAsFactors = FALSE)

############################
# File paths
############################

# Normalised methylation data (created in 01c_normalisation_RODAM.R)
beta_file  <- "Beta_Funnurm_NoSex_SNPs_NoCrosshybrRODAM_Adiponectin_Muhau_315.txt"
mval_file  <- "Mval_Funnurm_NoSex_SNPs_NoCrosshybrRODAM_Adiponectin_Muhau_n315.txt"
annot_file <- "Annot_NoXY_NoSNP_NoCrosshyb_Adiponectin_Muhau_n315.txt"

# Phenotype file with cell proportions
pheno_file <- "2024.10.31_Phenotype_RODAM_diabetes_case_n112_MM.csv"

############################
# Load data
############################

beta <- read.table(beta_file, sep = "\t", header = TRUE, check.names = FALSE)
mval <- read.table(mval_file, sep = "\t", header = TRUE, check.names = FALSE)
annot <- read.table(annot_file, sep = "\t", header = TRUE)

targets <- read.csv(pheno_file)
head(targets)

############################
# Subset: T2D cases only
############################

targets <- targets[targets$DM_Dichot == 1, ]
stopifnot(nrow(targets) == 112)

############################
# Define variables
############################

# Outcome
log_adip <- targets$log_Adiponectin

# Phenotypic covariates
sex  <- factor(targets$Sex)
age  <- targets$Age
BMI  <- targets$BMI
site <- factor(targets$Site)

# Technical covariates
batchArray <- factor(targets$BATCH_arry)
pos        <- targets$Sentrix_Position

# Cell proportions
CD8T <- targets$CD8T
CD4T <- targets$CD4T
NK   <- targets$NK
Bcell<- targets$Bcell
Mono <- targets$Mono
Gran <- targets$Gran

############################
# Match methylation data to samples
############################

ids  <- targets$Basename
ids2 <- paste0("X", ids)

mval_ft <- mval[, colnames(mval) %in% ids2]
beta_ft <- beta[, colnames(beta) %in% ids2]

############################
# Design matrices
############################

design_full <- model.matrix(
  ~ log_adip + sex + age +
    CD8T + CD4T + NK + Bcell + Mono + Gran +
    batchArray + pos + site + BMI
)

design_base <- model.matrix(~ log_adip)

############################
# EWAS: M-values (primary)
############################

fit_m <- lmFit(mval_ft, design_full)
fit_m <- eBayes(fit_m)

fit_m0 <- lmFit(mval_ft, design_base)
fit_m0 <- eBayes(fit_m0)

top_m_fdr <- topTable(
  fit_m,
  coef = "log_adip",
  number = nrow(mval_ft),
  adjust = "fdr",
  confint = TRUE
)

top_m_raw <- topTable(
  fit_m0,
  coef = "log_adip",
  number = nrow(mval_ft),
  adjust = "none",
  confint = TRUE
)

top_m_fdr <- top_m_fdr[order(rownames(top_m_fdr)), ]
top_m_raw <- top_m_raw[order(rownames(top_m_raw)), ]

############################
# Add annotation & SE (M-values)
############################

annot_match <- annot[match(rownames(top_m_fdr), rownames(annot)), ]
delta_mval  <- top_m_raw$logFC
SE_mval     <- (top_m_fdr$CI.R - top_m_fdr$CI.L) / 3.92

top_m_out <- cbind(
  top_m_fdr,
  DeltaMval = delta_mval,
  SE_Mval   = SE_mval
)

############################
# EWAS: Beta values (secondary)
############################

fit_b <- lmFit(beta_ft, design_full)
fit_b <- eBayes(fit_b)

fit_b0 <- lmFit(beta_ft, design_base)
fit_b0 <- eBayes(fit_b0)

top_b_fdr <- topTable(
  fit_b,
  coef = "log_adip",
  number = nrow(beta_ft),
  adjust = "fdr",
  confint = TRUE
)

top_b_raw <- topTable(
  fit_b0,
  coef = "log_adip",
  number = nrow(beta_ft),
  adjust = "none",
  confint = TRUE
)

top_b_fdr <- top_b_fdr[order(rownames(top_b_fdr)), ]
top_b_raw <- top_b_raw[order(rownames(top_b_raw)), ]

delta_beta <- top_b_raw$logFC

top_b_out <- cbind(
  top_b_fdr,
  DeltaBeta = delta_beta
)

############################
# Combine M-values and Betas
############################

top_b_out <- top_b_out[match(rownames(top_m_out), rownames(top_b_out)), ]

final_table <- cbind(
  top_m_out,
  top_b_out,
  annot_match
)

############################
# Save output
############################

write.table(
  final_table,
  "2025.03.16_TopTable_RODAM_Adiponectin_T2Dcases_n112.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = TRUE
)

############################
# Cleanup
############################

rm(list = ls())
gc()
