############################################
# Analysis: clusterProfiler enrichment
# Cohort: AADM
# Author: Muhulo Muhau Mungamba
############################################

rm(list = ls())
gc()

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(stringr)

aadm <- read.delim(
  "TopTables/2024.10.22_TopTable_AADM_Adiponectin_n593_age_sex_cells_bmi_t2d.txt",
  sep = "\t",
  header = TRUE
)

genes <- aadm %>%
  filter(P.Value < 1e-3 & !is.na(UCSC_RefGene_Name)) %>%
  pull(UCSC_RefGene_Name) %>%
  str_split(";") %>% unlist() %>% unique()

gene_df <- bitr(genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)

go <- enrichGO(gene_df$ENTREZID, org.Hs.eg.db, ont = "BP", readable = TRUE)
kegg <- enrichKEGG(gene_df$ENTREZID, organism = "hsa")

write.table(as.data.frame(go),   "Results/AADM_clusterProfiler_GO.txt",   sep="\t", quote=FALSE)
write.table(as.data.frame(kegg), "Results/AADM_clusterProfiler_KEGG.txt", sep="\t", quote=FALSE)

############################
# END OF SCRIPT
############################
