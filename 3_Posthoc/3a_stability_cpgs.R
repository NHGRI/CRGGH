##################
##Stability CpGs##
##################

# Load libraries
library(ggplot2)
library(ggbreak)
library(here)
library(knitr)
library(tidyverse)
library(tableone)
library(dplyr)
library(ggpubr)


# Getting the data

## Betas and Mvals extracted for CpGs at baseline and follow-up

rm(list=ls()) #Remove any existing objects in R 

betaBL <- read_tsv(here("Output/Extracted_CpGs", "Beta_per_sample_baseline.txt"))
betaBL <- as_tibble(betaBL)
nrow(betaBL)

betaFU <- read_tsv(here("Output/Extracted_CpGs", "Beta_per_sample_followup.txt"))
betaFU <- as_tibble(betaFU)
nrow(betaFU)


## Harmonize IDs

# load array info received from Eva
array <- read_csv(here("Data", "ArrayInformation.csv"))
colnames(array)

myvars <- c("Basename", "RodamID")
arrayx <- array[myvars]
colnames(arrayx)

# Create a lookup table from arrayx
lookup <- arrayx %>%
  mutate(Basename = paste0("X", Basename))


# Get the current column names of betaBL
current_colnames <- colnames(betaBL)

# Replace the matching Basename with RodamID
new_colnames <- current_colnames %>%
  map_chr(~ {
    match <- lookup$RodamID[lookup$Basename == .x]
    if(length(match) == 0) {
      .x  # Keep original name if no match is found
    } else {
      match
    }
  })

# Assign the new column names to betaBL
colnames(betaBL) <- new_colnames

## Transpose data frame BL Beta
betaBL_transposed <- t(betaBL)
betaBL_transposed_df <- as.data.frame(betaBL_transposed)
colnames(betaBL_transposed_df) <- betaBL_transposed_df[1, ]  # Set first row as column names
betaBL_transposed_df <- betaBL_transposed_df[-1, ]  # Remove the first row
betaBL_transposed_df <- tibble::rownames_to_column(betaBL_transposed_df, var = "RodamID")
betaBL <- betaBL_transposed_df


# Get the current column names of betaFU
current_colnames <- colnames(betaFU)

# Replace the matching Basename with cohort ID
new_colnames <- current_colnames %>%
  map_chr(~ {
    match <- lookup$RodamID[lookup$Basename == .x]
    if(length(match) == 0) {
      .x  # Keep original name if no match is found
    } else {
      match
    }
  })

# Assign the new column names to betaFU
colnames(betaFU) <- new_colnames

## Transpose data frame follow-up Beta
betaFU_transposed <- t(betaFU)
betaFU_transposed_df <- as.data.frame(betaFU_transposed)
colnames(betaFU_transposed_df) <- betaFU_transposed_df[1, ]  # Set first row as column names
betaFU_transposed_df <- betaFU_transposed_df[-1, ]  # Remove the first row
betaFU_transposed_df <- tibble::rownames_to_column(betaFU_transposed_df, var = "RodamID")
betaFU <- betaFU_transposed_df


## Merge baseline and follow-up data
merged_Beta <- merge(betaFU, betaBL, by = "ID", suffixes = c(".FU", ".BL"))
ncol(merged_Beta)
nrow(merged_Beta)

merged_Beta <- merged_Beta %>%
  mutate(across(-1, ~ as.numeric(as.character(.))))



# Plot stability over time

# Reshape data from wide to long format and convert to percentages
long_data <- merged_Beta %>%
  pivot_longer(cols = -RodamID, 
               names_to = c("CpG_site", "Timepoint"),
               names_sep = "\\.",
               values_to = "Proportion") %>%
  mutate(Proportion = Proportion * 100)  # Convert proportions to percentages

# Calculate the mean methylation percentage for each CpG site at baseline and follow-up
mean_data <- long_data %>%
  group_by(CpG_site, Timepoint) %>%
  summarise(mean_percentage = mean(Proportion, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Timepoint, values_from = mean_percentage) %>%
  drop_na()  


## Bland Altman plot
mean_data <- mean_data %>%
  mutate(average = (BL + FU) / 2,
         difference = FU - BL)

BAplot <- ggplot(mean_data, aes(x = average, y = difference)) +
  geom_point(color = "black") +  # Default points
  geom_point(data = subset(mean_data, CpG_site == "cg00036588"), 
             aes(x = average, y = difference), 
             color = "firebrick", size = 3) +  # Highlight cg00036588
  geom_point(data = subset(mean_data, CpG_site == "cg16759041"), 
             aes(x = average, y = difference), 
             color = "skyblue3", size = 3) +  # Highlight cg16759041
  geom_hline(yintercept = 0, linetype = "dashed") +  # Line at y=0
  labs(x = "% Mean Methylation (Timepoint 1 & 2)", 
       y = "% Methylation Difference") +
  theme_minimal() +
  theme(plot.background = element_rect(fill = "white", color = NA),  # White background
        axis.text.x = element_text(size = 16),  # Increase x-axis tick size
        axis.text.y = element_text(size = 16),  # Increase y-axis tick size
        axis.title.x = element_text(size = 18, face = "bold"),  # Increase x-axis title size
        axis.title.y = element_text(size = 18, face = "bold"))  # Increase y-axis title size

BAplot


