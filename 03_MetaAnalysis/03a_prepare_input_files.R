############################################################
# Title: 03a_prepare_input_files.R
# Project: Adiponectin EWAS – RODAM & AADM
# Author: Muhulo Muhau Mungamba
# Date: 30/10/2024
#
# Purpose:
# Prepare harmonized METAL input files from stratified EWAS
# results.
#
# This script prepares input only – no meta-analysis is run here.
############################################################

rm(list = ls())
gc()

library(data.table)
library(dplyr)

############################################################
# USER-DEFINED PATHS
############################################################

# Folder with EWAS toptables
ewas_dir <- "/project/MUHAU_RODAM/EWAS_Adiponectin/Results/EWAS_Toptables/"

# Output folder for METAL input files
out_dir <- "/project/MUHAU_RODAM/EWAS_Adiponectin/03_MetaAnalysis/METAL_inputs/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
# HELPER FUNCTION
############################################################

prepare_metal_input <- function(
    file,
    sample_size,
    cohort,
    t2d_status,
    output_name
) {
  
  message("Processing: ", output_name)
  
  tt <- fread(file)
  
  required_cols <- c("CpG", "chr", "pos", "logFCMval", "SE_mval", "P.Value.Mval")
  missing_cols <- setdiff(required_cols, colnames(tt))
  if (length(missing_cols) > 0) {
    stop("Missing columns in ", file, ": ", paste(missing_cols, collapse = ", "))
  }
  
  metal_df <- tt %>%
    mutate(
      CHR    = as.numeric(gsub("chr", "", chr)),
      POS    = as.numeric(pos),
      MARKER = paste(CHR, POS, sep = ":"),
      CPG    = CpG,
      BETA   = logFCMval,
      SE     = SE_mval,
      PVALUE = P.Value.Mval,
      N      = sample_size,
      COHORT = cohort,
      T2D    = t2d_status
    ) %>%
    select(MARKER, CPG, CHR, POS, BETA, SE, PVALUE, N)
  
  fwrite(
    metal_df,
    file = file.path(out_dir, output_name),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
}

############################################################
# SAMPLE SIZES (FROM MANUSCRIPT)
############################################################

# RODAM
n_rodam_cases    <- 112
n_rodam_controls <- 203
n_rodam_total    <- 315

# AADM
n_aadm_cases     <- 277
n_aadm_controls  <- 316
n_aadm_total     <- 593

############################################################
# 1. RODAM – STRATIFIED EWAS
############################################################

prepare_metal_input(
  file        = paste0(ewas_dir, "TopTable_RODAM_Adiponectin_T2Dcases.txt"),
  sample_size = n_rodam_cases,
  cohort      = "RODAM",
  t2d_status  = "Cases",
  output_name = "RODAM_T2Dcases_METAL.txt"
)

prepare_metal_input(
  file        = paste0(ewas_dir, "TopTable_RODAM_Adiponectin_T2Dcontrols.txt"),
  sample_size = n_rodam_controls,
  cohort      = "RODAM",
  t2d_status  = "Controls",
  output_name = "RODAM_T2Dcontrols_METAL.txt"
)

############################################################
# 2. AADM – STRATIFIED EWAS
############################################################

prepare_metal_input(
  file        = paste0(ewas_dir, "TopTable_AADM_Adiponectin_T2Dcases.txt"),
  sample_size = n_aadm_cases,
  cohort      = "AADM",
  t2d_status  = "Cases",
  output_name = "AADM_T2Dcases_METAL.txt"
)

prepare_metal_input(
  file        = paste0(ewas_dir, "TopTable_AADM_Adiponectin_T2Dcontrols.txt"),
  sample_size = n_aadm_controls,
  cohort      = "AADM",
  t2d_status  = "Controls",
  output_name = "AADM_T2Dcontrols_METAL.txt"
)

############################################################
# 3. WITHIN-COHORT META-ANALYSIS INPUTS
############################################################

# RODAM cases + controls
prepare_metal_input(
  file        = paste0(ewas_dir, "TopTable_RODAM_Adiponectin_all.txt"),
  sample_size = n_rodam_total,
  cohort      = "RODAM",
  t2d_status  = "All",
  output_name = "RODAM_cases_controls_METAL.txt"
)

# AADM cases + controls
prepare_metal_input(
  file        = paste0(ewas_dir, "TopTable_AADM_Adiponectin_all.txt"),
  sample_size = n_aadm_total,
  cohort      = "AADM",
  t2d_status  = "All",
  output_name = "AADM_cases_controls_METAL.txt"
)

############################################################
# DONE
############################################################

message("✔ All METAL input files successfully created in:")
message(out_dir)

save.image(file = file.path(out_dir, "03a_prepare_input_files.RData"))
savehistory(file = file.path(out_dir, "03a_prepare_input_files.Rhistory"))

###########################################################
#End of analysis
##########################################################