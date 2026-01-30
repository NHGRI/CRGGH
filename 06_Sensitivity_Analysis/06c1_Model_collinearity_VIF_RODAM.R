############################################################
# Title: Model collinearity assessment – RODAM
# Project: Adiponectin EWAS – Sensitivity Analysis
#
# Author: Muhulo Muhau Mungamba
# Date: 2025-10-17
#
# Description:
# This script evaluates multicollinearity among covariates
# included in the RODAM EWAS model using:
#  - Pearson correlation (BMI vs log-adiponectin)
#  - Variance Inflation Factors (VIF)
############################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

library(dplyr)
library(car)
library(ggplot2)

# ---- Set working directory ----
setwd("/home/P083854/lkg/Rodam/MUHAU_RODAM/EWAS_Adiponectin/Manuscript/Publication_ready/eBioMedicine/Revision/")

# ---- Load phenotype file ----
phenofile <- "/home/P083854/lkg/Rodam/Karlijn/EWAS_Adiponectin/Pheno/20230105_RODAM_Adiponectin_Phenofile_wCells.csv"
rodam <- read.csv(phenofile)

# ---- Prepare variables ----
rodam <- rodam %>%
  mutate(
    log_adip = log_Adiponectin,
    sex = factor(Sex),
    site = factor(Site),
    batchArray = factor(BATCH_arry),
    pos = factor(Sentrix_Position)
  )

# ---- Correlation: BMI vs log(Adiponectin) ----
cor_test <- cor.test(rodam$log_adip, rodam$BMI, use = "complete.obs")
print(cor_test)

p_cor <- ggplot(rodam, aes(x = BMI, y = log_adip)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(
    x = "BMI",
    y = "log(Adiponectin)",
    title = "RODAM: BMI vs log(Adiponectin)"
  ) +
  theme_minimal()

# ---- Save plot ----
outdir <- "06_Sensitivity_Analysis/06c_Model_Collinearity/outputs"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(outdir, "Correlation_RODAM_BMI_Adiponectin.png"),
  plot = p_cor,
  width = 6,
  height = 4,
  dpi = 300
)

# ---- VIF model ----
vif_model <- lm(
  log_adip ~ Age + sex + BMI + site + batchArray + pos +
    CD8T + CD4T + NK + Bcell + Mono + Gran,
  data = rodam
)

vif_values <- vif(vif_model)

vif_df <- as.data.frame(vif_values)
vif_df$Covariate <- rownames(vif_df)
rownames(vif_df) <- NULL

colnames(vif_df) <- c("GVIF", "Df", "GVIF_adj", "Covariate")

vif_df <- vif_df %>%
  mutate(
    VIF_adj = round(GVIF_adj, 2),
    Interpretation = case_when(
      VIF_adj < 5 ~ "No collinearity",
      VIF_adj < 10 ~ "Moderate collinearity",
      TRUE ~ "High collinearity"
    )
  ) %>%
  select(Covariate, VIF_adj, Interpretation)

print(vif_df)

write.csv(
  vif_df,
  file.path(outdir, "VIF_RODAM.csv"),
  row.names = FALSE
)

############################################################
# END OF SCRIPT
############################################################
