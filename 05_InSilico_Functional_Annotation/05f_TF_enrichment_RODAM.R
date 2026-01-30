############################################
# Analysis: Transcription factor enrichment
# Cohort: RODAM
#Author: Muhulo Muhau Mungamba
############################################

rm(list = ls())
gc()

library(enrichR)
library(dplyr)
library(stringr)

rodam <- read.delim(
  "TopTables/2025.10.14_TopTable_RODAM_Adiponectin_combined_T2Dcasescontrols_n315_age_sex_cells_pos_batcharry_BMI_site_DM.txt",
  sep = "\t",
  header = TRUE
)

genes <- rodam %>%
  filter(P.Value.Mval < 1e-4 & !is.na(UCSC_RefGene_Name)) %>%
  pull(UCSC_RefGene_Name) %>%
  str_split(";") %>% unlist() %>% unique()

tf <- enrichr(
  genes,
  databases = c("ChEA_2016", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
)

write.table(tf[["ChEA_2016"]], "Results/RODAM_TF_ChEA.txt", sep="\t", quote=FALSE)
write.table(tf[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]],
            "Results/RODAM_TF_ENCODE.txt", sep="\t", quote=FALSE)

############################
# END OF SCRIPT
############################
