############################################################
# RODAM Adiponectin EWAS
# Step 1: Phenotype preparation
#
# Purpose:
# - Prepare EWAS phenotype file for adiponectin analysis
# - Restrict to participants with adiponectin and DNAm data
# - Remove extreme adiponectin outliers on raw scale
# - Log-transform adiponectin
#
# Final output:
# - EWAS phenotype file (n = 315)
#
# Author: Muhau M. Mungamba
# Date: 10/10/2024
############################################################

## ---------------------------------------------------------
## 0. Housekeeping
## ---------------------------------------------------------

rm(list = ls())
gc()

## ---------------------------------------------------------
## 1. Load required packages
## ---------------------------------------------------------

library(dplyr)
library(haven)

## ---------------------------------------------------------
## 2. Read phenotype data
## ---------------------------------------------------------

# Path to prepared RODAM phenotype dataset
pheno_path <- "path/to/Rodam_Complete_Dataset_updateMay2022_withAdipokines_withEpigenID_JW.sav"

rodam <- read_sav(pheno_path)

## ---------------------------------------------------------
## 3. Restrict to participants with required data
## ---------------------------------------------------------

# Keep participants with adiponectin measurements
rodam_adip <- rodam %>%
  filter(!is.na(Adiponectin))

# Keep participants with DNA methylation data available
rodam_adip_dnAm <- rodam_adip %>%
  filter(!is.na(Epigenetics))

## ---------------------------------------------------------
## 4. Remove extreme adiponectin outliers (raw scale, ±3 SD)
## ---------------------------------------------------------

mean_adip <- mean(rodam_adip_dnAm$Adiponectin, na.rm = TRUE)
sd_adip   <- sd(rodam_adip_dnAm$Adiponectin, na.rm = TRUE)

rodam_clean <- rodam_adip_dnAm %>%
  filter(
    Adiponectin >= (mean_adip - 3 * sd_adip),
    Adiponectin <= (mean_adip + 3 * sd_adip)
  )

## ---------------------------------------------------------
## 5. Log-transform adiponectin
## ---------------------------------------------------------

rodam_clean <- rodam_clean %>%
  mutate(
    log_Adiponectin = log(Adiponectin)
  )

## Check final sample size
nrow(rodam_clean)  # Expected: 315

## ---------------------------------------------------------
## 6. Select EWAS-relevant variables
## ---------------------------------------------------------

rodam_ewas_pheno <- rodam_clean %>%
  select(
    RodamID,
    Sentrix_ID,
    Sentrix_Position,
    BATCH_arry,
    Batch_BSR,
    Sex,
    Age,
    BMI,
    Site,
    DM_Dichot,
    Adiponectin,
    log_Adiponectin
  )

## ---------------------------------------------------------
## 7. Save EWAS phenotype file
## ---------------------------------------------------------

write.csv(
  rodam_ewas_pheno,
  file = "RODAM_Adiponectin_EWAS_phenofile_n315.csv",
  row.names = FALSE,
  quote = FALSE
)

## ---------------------------------------------------------
## End of script
## ---------------------------------------------------------
