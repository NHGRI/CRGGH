############################################
# Title: QQ plots and lambda for Meta-analysis EWAS
# Project: Adiponectin EWAS (RODAM & AADM)
# Author: Muhulo Muhau Mungamba
# Date: 23/03/2025
#
# Description:
# This script generates QQ plots and calculates genomic inflation
# factors (lambda) for meta-analysis EWAS results using both raw
# and BACON-corrected p-values.
#
# Analyses included:
# 1. West Africans (RODAM + AADM) – cases + controls (MAIN)
# 2. West Africans – T2D cases
# 3. West Africans – T2D controls
#
# Notes:
# - Meta-analysis performed using METAL (Linux command line)
# - This script is for downstream QC and visualization in R
############################################

rm(list = ls())
gc()

############################
## Function: QQ plot
############################
qqplot_ewas <- function(pvector, ...) {
  if (!is.numeric(pvector)) stop("Input must be numeric.")
  
  pvector <- pvector[
    !is.na(pvector) &
      is.finite(pvector) &
      pvector > 0 &
      pvector < 1
  ]
  
  o <- -log10(sort(pvector))
  e <- -log10(ppoints(length(pvector)))
  
  plot(
    e, o,
    pch = 20,
    xlab = expression(Expected~~-log[10](italic(p))),
    ylab = expression(Observed~~-log[10](italic(p))),
    ...
  )
  
  abline(0, 1, col = "#B41A21", lwd = 3)
}

############################
## Function: lambda
############################
calc_lambda <- function(pvector) {
  pvector <- pvector[
    !is.na(pvector) &
      is.finite(pvector) &
      pvector > 0 &
      pvector < 1
  ]
  chisq <- qchisq(1 - pvector, 1)
  median(chisq) / qchisq(0.5, 1)
}

############################
## 1. West Africans (cases + controls)
############################

meta_WA_all <- read.csv(
  "Adiponectin_RODAM_AADM_T2D_cases_controls_metaanalysisN1_wFDR.txt"
)

# QQ plot – raw p-values
tiff("QQplot_Meta_WestAfricans_All_raw.tiff", width = 1500, height = 1500)
par(mai = c(1.6, 2.5, 0.5, 0.5), mgp = c(7, 2, 0))
qqplot_ewas(meta_WA_all$P.value, ylim = c(0, 12))
dev.off()

# QQ plot – Bacon-corrected
tiff("QQplot_Meta_WestAfricans_All_BACON.tiff", width = 1500, height = 1500)
par(mai = c(1.6, 2.5, 0.5, 0.5), mgp = c(7, 2, 0))
qqplot_ewas(meta_WA_all$BaconPvalue, ylim = c(0, 12))
dev.off()

lambda_WA_raw   <- calc_lambda(meta_WA_all$P.value)
lambda_WA_bacon <- calc_lambda(meta_WA_all$BaconPvalue)

print(paste("Lambda WA all (raw):", round(lambda_WA_raw, 3)))
print(paste("Lambda WA all (BACON):", round(lambda_WA_bacon, 3)))

############################
## 2. West Africans – T2D cases
############################

meta_WA_cases <- read.csv(
  "Adiponectin_RODAM_AADM_T2D_cases_metaanalysisN1_wFDR.txt"
)

tiff("QQplot_Meta_WestAfricans_T2Dcases_BACON.tiff", width = 1500, height = 1500)
par(mai = c(1.6, 2.5, 0.5, 0.5), mgp = c(7, 2, 0))
qqplot_ewas(meta_WA_cases$BaconPvalue, ylim = c(0, 12))
dev.off()

lambda_WA_cases <- calc_lambda(meta_WA_cases$BaconPvalue)
print(paste("Lambda WA T2D cases (BACON):", round(lambda_WA_cases, 3)))

############################
## 3. West Africans – T2D controls
############################

meta_WA_controls <- read.csv(
  "Adiponectin_RODAM_AADM_T2D_controls_metaanalysisN1_wFDR.txt"
)

tiff("QQplot_Meta_WestAfricans_T2Dcontrols_BACON.tiff", width = 1500, height = 1500)
par(mai = c(1.6, 2.5, 0.5, 0.5), mgp = c(7, 2, 0))
qqplot_ewas(meta_WA_controls$BaconPvalue, ylim = c(0, 12))
dev.off()

lambda_WA_controls <- calc_lambda(meta_WA_controls$BaconPvalue)
print(paste("Lambda WA T2D controls (BACON):", round(lambda_WA_controls, 3)))

############################
## Save workspace
############################

save.image("03g_QQplots_MetaAnalysis.RData")

############################
## END OF SCRIPT
############################
