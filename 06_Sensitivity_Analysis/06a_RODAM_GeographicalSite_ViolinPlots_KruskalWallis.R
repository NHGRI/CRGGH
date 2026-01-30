############################################################
# Title: Geographical site-specific variation in DNA methylation
# Project: Adiponectin EWAS – Sensitivity Analysis
# Figure: Figure 6 (Violin plots by site)
#
# Author: Muhulo Muhau Mungamba
# Date: 2025-06-16
#
# Description:
# Sensitivity analysis assessing geographical variation in DNA
# methylation levels at genome-wide significant CpG sites across
# RODAM recruitment sites (Rural Ghana, Urban Ghana, Amsterdam).
#
# Statistical approach:
# - Visualization using violin plots
# - Global differences tested using Kruskal–Wallis test
# - Post-hoc pairwise Wilcoxon tests with FDR correction
#
# Notes:
# - Analysis restricted to RODAM (site information not available in AADM)
# - Beta-values used for visualization
############################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

# ---- Load libraries ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(tibble)

# ---- Set working directory (EDIT AS NEEDED) ----
setwd("/project/MUHAU_RODAM/")

# ---- Load phenotype data ----
phenadip <- read.csv(
  "/project/Karlijn/EWAS_Adiponectin/Pheno/20230105_RODAM_Adiponectin_Phenofile_wCells.csv"
)

# Ensure Basename matches methylation matrices
phenadip$Basename <- paste0("X", phenadip$Basename)

# ---- Load methylation data (significant CpGs only) ----
betas <- read.table("betas_adiponectin_significantCpGs_MM.txt", header = TRUE)
mvals <- read.table("Mvals_adiponectin_significantCpGs_MM.txt", header = TRUE)

# ---- Reshape beta-values ----
betas_long <- betas %>%
  pivot_longer(cols = -CPG, names_to = "Basename", values_to = "Beta") %>%
  pivot_wider(names_from = CPG, values_from = Beta)

# ---- Reshape M-values ----
mvals_df <- as.data.frame(t(mvals))
mvals_df <- rownames_to_column(mvals_df, var = "Basename")

# ---- Merge phenotype + methylation ----
df <- phenadip %>%
  left_join(betas_long, by = "Basename") %>%
  left_join(mvals_df, by = "Basename")

# ---- Convert beta-values to percentages ----
df <- df %>%
  mutate(
    cg03546163_bp = cg03546163.x * 100,
    cg02561343_bp = cg02561343.x * 100,
    cg23969380_bp = cg23969380.x * 100
  )

# ---- Recode site labels ----
df$Site <- factor(
  df$Site,
  levels = c(1, 4, 5),
  labels = c("Rural Ghana", "Urban Ghana", "Amsterdam")
)

table(df$Site)

# ==========================================================
# FUNCTION: Violin plot + Kruskal–Wallis test
# ==========================================================
plot_violin_kw <- function(data, outcome, gene, cpg) {
  
  kw <- kruskal.test(as.formula(paste(outcome, "~ Site")), data = data)
  pval <- signif(kw$p.value, 3)
  
  ggplot(data, aes(x = Site, y = .data[[outcome]], fill = Site)) +
    geom_violin(trim = FALSE) +
    stat_summary(fun = median, geom = "point", size = 3, color = "black") +
    labs(
      title = bquote(italic(.(gene)) ~ "(" * .(cpg) * ")"),
      x = "Site",
      y = "Percentage methylation"
    ) +
    annotate(
      "text",
      x = 3.3,
      y = max(data[[outcome]], na.rm = TRUE) * 0.98,
      label = paste0("Kruskal–Wallis p = ", pval),
      size = 5
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 15),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

# ---- Generate plots ----
plot_violin_kw(df, "cg03546163_bp", "FKBP5", "cg03546163")
plot_violin_kw(df, "cg02561343_bp", "UST", "cg02561343")
plot_violin_kw(df, "cg23969380_bp", "ADGRD1", "cg23969380")

# ==========================================================
# Pairwise Wilcoxon post-hoc tests (FDR corrected)
# ==========================================================
run_pairwise <- function(data, outcome, gene) {
  pw <- pairwise.wilcox.test(
    data[[outcome]],
    data$Site,
    p.adjust.method = "fdr"
  )
  
  out <- as.data.frame(as.table(pw$p.value))
  colnames(out) <- c("Group1", "Group2", "Adjusted_P_Value")
  out$CpG <- outcome
  out$Gene <- gene
  out
}

pairwise_results <- bind_rows(
  run_pairwise(df, "cg03546163_bp", "FKBP5"),
  run_pairwise(df, "cg02561343_bp", "UST"),
  run_pairwise(df, "cg23969380_bp", "ADGRD1")
)

print(pairwise_results)

# Optional save
# write.csv(pairwise_results,
#           "RODAM_site_pairwise_Wilcoxon_results.csv",
#           row.names = FALSE)

############################################################
# END OF SCRIPT
############################################################
