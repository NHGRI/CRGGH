####################
##Temporality CpGs##
####################

# Load libraries
library(ggplot2)
library(ggbreak)
library(here)
library(tidyverse)
library(dplyr)


## Load data and merge
rm(list=ls())
dta <- read.table(here("Data", "Longitudinal_data.txt"), header = TRUE)
colnames(dta)
nrow(dta)

beta <- read.table(here("Output/Extracted_CpGs", "Beta_per_sample.txt"), header=TRUE)
nrow(beta)
colnames(beta)

dfm <- merge(dta, beta, by = "Basename")
colnames(dfm)


# Create tertiles for HbA1c
dfm <- dfm %>%
  mutate(HbA1c_tertile = ntile(HBA1C_2, 3))

#mean hbA1c per tertile
dfm %>%
  group_by(HbA1c_tertile) %>%
  summarise(mean_HbA1c = mean(HBA1C_2, na.rm = TRUE))

dfm %>%
  group_by(HbA1c_tertile) %>%
  summarise(median_beta = median(cpg, na.rm = TRUE))

kruskal.test(cpg ~ factor(HbA1c_tertile), data = dfm)

# Box plot for baseline methylation vs. HbA1c tertiles
plot2 <- ggplot(df2_nomed, aes(x = factor(HbA1c_tertile), y = cpg, fill = factor(HbA1c_tertile))) +
  geom_boxplot() +
  scale_fill_manual(values = c("firebrick", "skyblue3", "darkgreen"), labels = c("Tertile 1", "Tertile 2", "Tertile 3")) +  # Custom colors
  scale_x_discrete(labels = c("1" = "Tertile 1", "2" = "Tertile 2", "3" = "Tertile 3")) +  # Set x-axis labels
  labs(x = "HbA1c at Follow-up in Tertiles",
       y = "% Methylation cpg at Baseline") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 21, face = "bold"),
    axis.title.y = element_text(size = 21, face = "bold"),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.position = "none"
  )

print(plot2)


