############################################
# RODAM EWAS Adiponectin
# Manhattan plot
#
# Author: Muhulo Muhau Mungamba
# Date: 2025-03-25
# Analysis: RODAM combined (T2D cases + controls, n = 315)
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
# Load EWAS results (combined)
############################

toptables <- read.table(
  "2025.10.14_TopTable_RODAM_Adiponectin_combined_T2Dcasescontrols_n315_age_sex_cells_pos_batcharry_BMI_site.txt",
  header = TRUE,
  sep = "\t"
)

head(toptables)

############################
# Prepare data for Manhattan plot
############################

toptables <- toptables %>%
  mutate(
    chr = gsub("chr", "", chr),
    chr = as.numeric(chr),
    pos = as.numeric(pos),
    pval = BaconPvalue_Mvalues
  ) %>%
  filter(
    !is.na(chr),
    !is.na(pos),
    chr %in% 1:22,
    pos > 0,
    !is.na(pval),
    pval > 0,
    pval < 1
  )

############################
# Build Manhattan plot coordinates
############################

chrom_colors <- c("seagreen3", "skyblue3")

coor <- c()
cols <- c()
ticks <- c()
p <- c()

old_max <- 0

for (chr_i in 1:22) {
  chr_data <- toptables[toptables$chr == chr_i, ]
  
  if (nrow(chr_data) > 0) {
    chr_pos <- chr_data$pos + old_max
    
    coor <- c(coor, chr_pos)
    p <- c(p, chr_data$pval)
    cols <- c(cols, rep(chrom_colors[chr_i %% 2 + 1], nrow(chr_data)))
    ticks <- c(ticks, mean(chr_pos))
    
    old_max <- max(chr_pos)
  }
}

stopifnot(length(coor) == length(p))

############################
# Manhattan plot
############################

png(
  "2025.10.14_Manhattan_RODAM_combined_n315.png",
  width = 1600,
  height = 600,
  res = 120
)

plot(
  coor,
  -log10(p),
  pch = 20,
  col = cols,
  xaxt = "n",
  xlab = "Chromosome",
  ylab = expression(-log[10](italic(p))),
  cex = 0.6,
  cex.lab = 1.5,
  cex.axis = 1.3,
  ylim = c(0, max(-log10(p), na.rm = TRUE) + 1)
)

axis(1, at = ticks, labels = 1:22)

# Epigenome-wide significance threshold
abline(h = -log10(2e-7), col = "red", lty = 3, lwd = 2)

dev.off()

############################
# Clean up
############################

rm(toptables, coor, cols, ticks, p, old_max)
gc()
