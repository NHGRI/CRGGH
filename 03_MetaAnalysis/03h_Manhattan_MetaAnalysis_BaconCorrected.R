############################################
# Title: Manhattan plots for Meta-analysis EWAS (BACON corrected)
# Project: Adiponectin EWAS (RODAM & AADM)
# Author: Muhulo Muhau Mungamba
# Date: 23/03/2025
#
# Description:
# This script generates Manhattan plots using BACON-corrected
# p-values from METAL meta-analysis results.
#
# Analyses included:
# 1. West Africans (RODAM + AADM) – cases + controls (MAIN)
# 2. West Africans – T2D cases
# 3. West Africans – T2D controls
#
# Notes:
# - Meta-analysis performed using METAL (Linux)
# - BACON correction applied prior to plotting
############################################

rm(list = ls())
gc()

############################
## Load libraries
############################
library(dplyr)

############################
## Function: Manhattan plot
############################
manhattan_ewas <- function(df,
                           chr_col = "Chr",
                           pos_col = "Pos",
                           p_col,
                           outfile,
                           genomewide_line,
                           suggestive_line = NULL,
                           ylim_max = 15) {
  
  df <- df %>%
    mutate(
      chr = as.numeric(gsub("chr", "", .data[[chr_col]])),
      pos = as.numeric(.data[[pos_col]])
    ) %>%
    filter(chr %in% 1:22, !is.na(pos), pos > 0)
  
  df <- df[order(df$chr, df$pos), ]
  
  df$logp <- -log10(df[[p_col]])
  
  chr_lengths <- df %>%
    group_by(chr) %>%
    summarise(max_pos = max(pos)) %>%
    mutate(cumlen = cumsum(max_pos) - max_pos)
  
  df <- df %>%
    left_join(chr_lengths, by = "chr") %>%
    mutate(pos_cum = pos + cumlen)
  
  axis_df <- chr_lengths %>%
    mutate(center = cumlen + max_pos / 2)
  
  cols <- rep(c("seagreen3", "skyblue3"), length.out = 22)
  point_cols <- cols[df$chr]
  
  png(outfile, height = 600, width = 1600, res = 96)
  
  plot(
    df$pos_cum,
    df$logp,
    col = point_cols,
    pch = 20,
    xaxt = "n",
    xlab = "Chromosome",
    ylab = expression(-log[10](italic(p))),
    ylim = c(0, ylim_max),
    cex.lab = 1.5,
    cex.axis = 1.3
  )
  
  axis(1, at = axis_df$center, labels = axis_df$chr)
  
  abline(h = -log10(genomewide_line), col = "red", lty = 3)
  
  if (!is.null(suggestive_line)) {
    abline(h = -log10(suggestive_line), col = "blue", lty = 3)
  }
  
  dev.off()
}

############################
## 1. West Africans – cases + controls (MAIN)
############################

meta_WA_all <- read.csv(
  "Adiponectin_RODAM_AADM_T2D_cases_controls_metaanalysisN1_wFDR.txt"
)

manhattan_ewas(
  df = meta_WA_all,
  p_col = "BaconPvalue",
  outfile = "Manhattan_Meta_WestAfricans_All_BACON.png",
  genomewide_line = 4.46e-07
)

############################
## 2. West Africans – T2D cases
############################

meta_WA_cases <- read.csv(
  "Adiponectin_RODAM_AADM_T2D_cases_metaanalysisN1_wFDR.txt"
)

manhattan_ewas(
  df = meta_WA_cases,
  p_col = "BaconPvalue",
  outfile = "Manhattan_Meta_WestAfricans_T2Dcases_BACON.png",
  genomewide_line = 1.00e-08
)

############################
## 3. West Africans – T2D controls
############################

meta_WA_controls <- read.csv(
  "Adiponectin_RODAM_AADM_T2D_controls_metaanalysisN1_wFDR.txt"
)

manhattan_ewas(
  df = meta_WA_controls,
  p_col = "BaconPvalue",
  outfile = "Manhattan_Meta_WestAfricans_T2Dcontrols_BACON.png",
  genomewide_line = 2.00e-07,
  suggestive_line = 7.17e-06
)

############################
## Save workspace
############################
save.image("03h_Manhattan_MetaAnalysis_BaconCorrected.RData")

############################
## END OF SCRIPT
############################
