############################################################
# RODAM Adiponectin EWAS
# Step 1c: Normalisation and probe filtering
#
# Purpose:
# - Functional normalisation of DNAm data (funnorm)
# - Remove sex chromosome probes
# - Remove probes with SNPs and cross-hybridisation
# - Generate beta values, M-values, and annotation files
#
# Input:
# - IDAT files
# - Phenotype file with cell counts
#
# Output:
# - Normalised beta values
# - Normalised M-values
# - Probe annotation file
# - QC plots
#
# Author: Muhulo Muhau Mungamba
# Date: 15/10/2024
############################################################

## ---------------------------------------------------------
## 0. Housekeeping
## ---------------------------------------------------------

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

## ---------------------------------------------------------
## 1. Load required packages
## ---------------------------------------------------------

library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(minfiData)
library(limma)
library(affy)

## ---------------------------------------------------------
## 2. Define paths
## ---------------------------------------------------------

# Working directory (adjust to your system)
workdirectory <- "/path/to/01_RODAM/normalisation/"
setwd(workdirectory)

# IDAT path
idatspath <- "/path/to/Rodam/idats/"

# Phenotype file with cell counts (from 01b)
pheno_file <- "RODAM_Adiponectin_EWAS_phenofile_n315_wCells.csv"

# Cross-hybridisation probe list (Chen et al. 2013; African-specific SNPs)
cross_file <- "/path/to/crosshybEnSNPSRODAM.txt"

## ---------------------------------------------------------
## 3. Read phenotype file
## ---------------------------------------------------------

targets <- read.csv(pheno_file, header = TRUE)
head(targets)

## ---------------------------------------------------------
## 4. Read IDAT files
## ---------------------------------------------------------

RGset <- read.metharray.exp(
  base = idatspath,
  targets = targets,
  force = TRUE
)

pd <- pData(RGset)

## ---------------------------------------------------------
## 5. Extract unnormalised beta values
## ---------------------------------------------------------

beta_raw <- getBeta(RGset)  # ~485k probes (450k array)

## ---------------------------------------------------------
## 6. Functional normalisation
## ---------------------------------------------------------

GMset <- preprocessFunnorm(
  RGset,
  nPCs = 2,
  sex = c("male", "female"),
  bgCorr = TRUE,
  dyeCorr = TRUE,
  verbose = TRUE
)

## ---------------------------------------------------------
## 7. Probe annotation and filtering
## ---------------------------------------------------------

annotation <- getAnnotation(GMset)
probe.features <- as.matrix(annotation)

## Remove probes on sex chromosomes
sex_idx <- which(annotation$chr %in% c("chrX", "chrY"))
GMset_noSex <- GMset[-sex_idx, ]

## Remove SNP-affected and cross-hybridising probes
probes4removal <- read.csv(
  cross_file,
  stringsAsFactors = FALSE,
  col.names = 1
)[, 1]

GMset_culled <- GMset_noSex[
  !featureNames(GMset_noSex) %in% probes4removal,
]

gc()

## ---------------------------------------------------------
## 8. Generate beta values and M-values
## ---------------------------------------------------------

mval <- getM(GMset_culled)
beta <- getBeta(GMset_culled)

## Remove probes with infinite M-values
inf_idx <- which(mval == -Inf, arr.ind = TRUE)
bad_probes <- unique(rownames(inf_idx))

mval_clean <- mval[!rownames(mval) %in% bad_probes, ]
beta_clean <- beta[!rownames(beta) %in% bad_probes, ]

gc()

## ---------------------------------------------------------
## 9. Save output files
## ---------------------------------------------------------

# Unnormalised beta values
write.table(
  beta_raw,
  "Betaraw_RODAM_Adiponectin_n315.txt",
  sep = "\t",
  quote = FALSE
)

# Normalised M-values
write.table(
  mval_clean,
  "Mval_Funnorm_NoSex_SNPs_NoCrosshybr_RODAM_Adiponectin_n315.txt",
  sep = "\t",
  quote = FALSE
)

# Normalised beta values
write.table(
  beta_clean,
  "Beta_Funnorm_NoSex_SNPs_NoCrosshybr_RODAM_Adiponectin_n315.txt",
  sep = "\t",
  quote = FALSE
)

# Annotation file (matching filtered probes)
annot_out <- data.frame(
  probe.features[match(rownames(beta_clean), rownames(probe.features)), ]
)

write.table(
  annot_out,
  "Annot_NoXY_NoSNP_NoCrosshybr_RODAM_Adiponectin_n315.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

## ---------------------------------------------------------
## 10. Quality control plots
## ---------------------------------------------------------

# Raw beta density
pdf("QC_densityplot_betaraw_RODAM_Adiponectin.pdf")
densityPlot(RGset, main = "Raw beta values", xlab = "Beta")
dev.off()

# Raw beta bean plot
pdf("QC_densitybean_betaraw_RODAM_Adiponectin.pdf")
densityBeanPlot(RGset, sampNames = pd$Sample_ID, sampGroups = pd$Site)
dev.off()

# Bisulfite conversion controls
pdf("QC_controlstrip_betaraw_RODAM_Adiponectin.pdf")
controlStripPlot(
  RGset,
  controls = "BISULFITE CONVERSION II",
  sampNames = pd$Sample_ID
)
dev.off()

# Normalised beta density
pdf("QC_densityplot_beta_funnorm_RODAM_Adiponectin.pdf")
densityPlot(beta_clean, main = "Normalised beta values", xlab = "Beta")
dev.off()

# Normalised beta bean plot
pdf("QC_densitybean_beta_funnorm_RODAM_Adiponectin.pdf")
densityBeanPlot(beta_clean, sampNames = pd$Sample_ID, sampGroups = pd$Site)
dev.off()

## ---------------------------------------------------------
## End of script
## ---------------------------------------------------------
