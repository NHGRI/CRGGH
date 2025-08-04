#############################
##Prepare Methylation Input##
#############################

# load libraries
library(dplyr)

# set working directory
setwd("/QTLtools")

# Read M-values
m_values <- read.table("/Data/Mvalues.txt", header = TRUE)
nrow(m_values)

# Read the annotation file
annotation <- read.table("/Data/Annotation.txt", header = TRUE, sep = "\t")
nrow(annotation)

# Subset necessary columns from the annotation file and set missing UCSC_RefGene_Name to "NA"
annotation_subset <- annotation %>%
  select(CpG, strand, chr, start, end, UCSC_RefGene_Name) %>%
  mutate(UCSC_RefGene_Name = ifelse(UCSC_RefGene_Name == "", "NA", UCSC_RefGene_Name))

# Merge annotation with M-values
bed_data <- merge(annotation_subset, m_values, by = "CpG", all.x=TRUE)
nrow(bed_data)

# Reorder columns to match the required BED format
bed_data_formatted <- bed_data %>%
  select(chr, start, end, CpG, UCSC_RefGene_Name, strand, everything()) %>%
  rename(pid = CpG, gid = UCSC_RefGene_Name)

# Write to a new file in BED format
write.table(bed_data_formatted, "Mvalues.bed", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
