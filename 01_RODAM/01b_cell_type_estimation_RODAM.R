############################################################
# RODAM Adiponectin EWAS
# Step 1b: Blood cell-type proportion estimation
#
# Purpose:
# - Estimate blood cell-type proportions using the Houseman method
# - Merge cell-type proportions with EWAS phenotype file
#
# Input:
# - RODAM_Adiponectin_EWAS_phenofile_n315.csv
#
# Output:
# - RODAM_Adiponectin_EWAS_phenofile_n315_wCells.csv
#
# Author: Muhulo Muhau Mungamba
# Based on original script by: Karlijn Meeks (2023)
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
library(FlowSorted.Blood.450k)

## ---------------------------------------------------------
## 2. Read phenotype file
## ---------------------------------------------------------

targets <- read.csv(
  "RODAM_Adiponectin_EWAS_phenofile_n315.csv",
  header = TRUE
)

head(targets)

## ---------------------------------------------------------
## 3. Read IDAT files
## ---------------------------------------------------------

# Path to IDAT files
idatspath <- "path/to/Rodam/idats"
BasePath <- idatspath

RGset <- read.metharray.exp(
  base = BasePath,
  targets = targets
)

pd <- pData(RGset)

## ---------------------------------------------------------
## 4. Estimate blood cell-type proportions
## ---------------------------------------------------------

cell_est <- estimateCellCounts(
  RGset,
  compositeCellType = "Blood",
  cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"),
  returnAll = TRUE,
  meanPlot = FALSE,
  verbose = FALSE
)

## ---------------------------------------------------------
## 5. Prepare cell proportion table
## ---------------------------------------------------------

cell_counts <- cell_est$counts
cell_df <- as.data.frame(cell_counts)

# Add Basename column
cell_df$Basename <- rownames(cell_df)

# Reorder columns
cell_df <- cell_df[, c("Basename", setdiff(names(cell_df), "Basename"))]

head(cell_df)

## ---------------------------------------------------------
## 6. Merge cell proportions with phenotype file
## ---------------------------------------------------------

pheno_wCells <- merge(
  targets,
  cell_df,
  by = "Basename",
  all.x = TRUE
)

head(pheno_wCells)

## ---------------------------------------------------------
## 7. Save updated phenotype file
## ---------------------------------------------------------

write.table(
  pheno_wCells,
  file = "RODAM_Adiponectin_EWAS_phenofile_n315_wCells.csv",
  sep = ",",
  quote = FALSE,
  row.names = FALSE
)

## ---------------------------------------------------------
## End of script
## ---------------------------------------------------------
