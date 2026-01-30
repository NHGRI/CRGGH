############################################################
# Title: Model collinearity assessment – AADM
# Project: Adiponectin EWAS – Sensitivity Analysis
#
# Author: Muhulo Muhau Mungamba
# Date: 2025-11-21
############################################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

library(dplyr)
library(car)
library(ggplot2)

# ---- Set working directory ----
setwd("C:/Users/meekska/OneDrive - National Institutes of Health/H Drive/Supervision/Muhau Mungamba/Adiponectin EWAS/EBiomedicine_Revision/")

# ---- Load phenotype file ----
phenofile <- "2024.10.22_adiponectin_model_n593_wCells.csv"
aadm <- read.csv(phenofile)

# ---- Prepare variables ----
aadm <- aadm %>%
  mutate(
    log_adip = log_adiponectin,
    sex = factor(sex)
  )

# ---- Correlation: BMI vs log(Adiponectin) ----
cor_test <- cor.test(aadm$log_adip, aadm$bmi, use = "complete.obs")
print(cor_test)

p_cor <- ggplot(aadm, aes(x = bmi, y = log_adip)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(
    x = "BMI",
    y = "log(Adiponectin)",
    title = "AADM: BMI vs log(Adiponectin)"
  ) +
  theme_minimal()

# ---- Save plot ----
outdir <- "06_Sensitivity_Analysis/06c_Model_Collinearity/outputs"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(outdir, "Correlation_AADM_BMI_Adiponectin.png"),
  plot = p_cor,
  width = 6,
  height = 4,
  dpi = 300
)

# ---- VIF model ----
vif_model <- lm(
  log_adip ~ age + sex + bmi + CD8T + CD4T + NK + Bcell + Mono + Neu,
  data = aadm
)

vif_values <- vif(vif_model)

vif_df <- data.frame(
  Covariate = names(vif_values),
  VIF = round(vif_values, 2)
)

vif_df$Interpretation <- ifelse(
  vif_df$VIF < 5, "No collinearity",
  ifelse(vif_df$VIF < 10, "Moderate collinearity", "High collinearity")
)

print(vif_df)

write.csv(
  vif_df,
  file.path(outdir, "VIF_AADM.csv"),
  row.names = FALSE
)

############################################################
# END OF SCRIPT
############################################################
