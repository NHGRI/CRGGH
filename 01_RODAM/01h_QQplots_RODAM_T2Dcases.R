############################################
# RODAM EWAS Adiponectin
# QQ plots and lambda calculation
#
# Author: Muhulo Muhau Mungamba
# Date: 2025-03-25
# Analysis: RODAM T2D cases only (n = 112)
############################################

rm(list = ls())
gc()

############################
# Load libraries
############################

library(dplyr)

options(stringsAsFactors = FALSE)

############################
# Set working directory
############################
# EDIT if needed
setwd("/project/MUHAU_RODAM/EWAS_Adiponectin/Results/RODAM")

############################
# Load EWAS results (T2D cases)
############################

toptablesadip <- read.table(
  "2025.03.16_TopTable_RODAM_Adiponectin_T2Dcases_n112_age_sex_cells_pos_batcharry_BMI_site.txt",
  header = TRUE,
  sep = "\t"
)

head(toptablesadip)

############################
# QQ plot function
############################

qqplot_ewas <- function(pvector, ...) {
  if (!is.numeric(pvector))
    stop("Input must be numeric.")
  
  pvector <- pvector[
    !is.na(pvector) &
      !is.nan(pvector) &
      is.finite(pvector) &
      pvector > 0 &
      pvector < 1
  ]
  
  o <- -log10(sort(pvector, decreasing = FALSE))
  e <- -log10(ppoints(length(pvector)))
  
  def_args <- list(
    pch = 20,
    xlim = c(0, max(e)),
    ylim = c(0, max(o)),
    xlab = expression(Expected ~ -log[10](italic(p))),
    ylab = expression(Observed ~ -log[10](italic(p)))
  )
  
  dotargs <- list(...)
  
  do.call(
    "plot",
    c(
      list(x = e, y = o),
      def_args[!names(def_args) %in% names(dotargs)],
      dotargs
    )
  )
  
  abline(0, 1, col = "#B41A21", lwd = 3)
}

############################
# QQ plot – raw M-value p-values
############################

tiff(
  "2025.03.16_QQplot_RODAM_T2Dcases_age_sex_cells_batch_pos_BMI_site.tiff",
  width = 1500,
  height = 1500,
  res = 300
)

par(mai = c(1.6, 2.5, 0.5, 0.5), mgp = c(7, 2, 0))

qqplot_ewas(
  toptablesadip$P.Value.Mval,
  frame.plot = FALSE,
  cex = 1.8,
  cex.lab = 2.2,
  cex.axis = 1.8,
  bg = "white",
  ylim = c(0, 12)
)

dev.off()

############################
# QQ plot – bacon-corrected p-values
############################

tiff(
  "2025.03.16_QQplot_RODAM_T2Dcases_age_sex_cells_batch_pos_BMI_site_bacon.tiff",
  width = 1500,
  height = 1500,
  res = 300
)

par(mai = c(1.6, 2.5, 0.5, 0.5), mgp = c(7, 2, 0))

qqplot_ewas(
  toptablesadip$BaconPvalue_Mvalues,
  frame.plot = FALSE,
  cex = 1.8,
  cex.lab = 2.2,
  cex.axis = 1.8,
  bg = "white",
  ylim = c(0, 12)
)

dev.off()

############################
# Lambda calculation
############################

# Raw p-values
chisq_raw <- qchisq(1 - toptablesadip$P.Value.Mval, 1)
lambda_raw <- median(chisq_raw, na.rm = TRUE) / qchisq(0.5, 1)

# Bacon-corrected p-values
chisq_bacon <- qchisq(1 - toptablesadip$BaconPvalue_Mvalues, 1)
lambda_bacon <- median(chisq_bacon, na.rm = TRUE) / qchisq(0.5, 1)

lambda_results <- data.frame(
  Cohort = "RODAM",
  Group = "T2D cases",
  N = 112,
  Lambda_raw = lambda_raw,
  Lambda_bacon = lambda_bacon
)

print(lambda_results)

write.table(
  lambda_results,
  "2025.03.16_Lambda_RODAM_T2Dcases.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################
# Clean environment
############################

rm(
  chisq_raw,
  chisq_bacon,
  lambda_raw,
  lambda_bacon,
  toptablesadip
)

gc()
