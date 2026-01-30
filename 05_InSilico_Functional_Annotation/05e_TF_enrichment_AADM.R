############################################
# Analysis: Transcription factor enrichment
# Cohort: AADM
#Author: Muhulo Muhau Mungamba
############################################

rm(list = ls())
gc()

library(enrichR)
library(dplyr)
library(stringr)

aadm <- read.delim(
  "TopTables/2024.10.22_TopTable_AADM_Adiponectin_n593_age_sex_cells_bmi_t2d.txt",
  sep = "\t",
  header = TRUE
)

genes <- aadm %>%
  filter(P.Value < 1e-4 & !is.na(UCSC_RefGene_Name)) %>%
  pull(UCSC_RefGene_Name) %>%
  str_split(";") %>% unlist() %>% unique()

tf <- enrichr(
  genes,
  databases = c("ChEA_2016", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
)

write.table(tf[["ChEA_2016"]], "Results/AADM_TF_ChEA.txt", sep="\t", quote=FALSE)
write.table(tf[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]],
            "Results/AADM_TF_ENCODE.txt", sep="\t", quote=FALSE)

############################
# END OF SCRIPT
############################
