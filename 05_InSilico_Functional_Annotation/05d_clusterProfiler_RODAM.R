############################################
# Analysis: clusterProfiler enrichment
# Cohort: RODAM
#Author: Muhulo Muhau Mungamba
############################################

rm(list = ls())
gc()

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(stringr)

rodam <- read.delim(
  "TopTables/2025.10.14_TopTable_RODAM_Adiponectin_combined_T2Dcasescontrols_n315_age_sex_cells_pos_batcharry_BMI_site_DM.txt",
  sep = "\t",
  header = TRUE
)

genes <- rodam %>%
  filter(P.Value.Mval < 1e-3 & !is.na(UCSC_RefGene_Name)) %>%
  pull(UCSC_RefGene_Name) %>%
  str_split(";") %>% unlist() %>% unique()

gene_df <- bitr(genes, "SYMBOL", "ENTREZID", org.Hs.eg.db)

go <- enrichGO(gene_df$ENTREZID, org.Hs.eg.db, ont = "BP", readable = TRUE)
kegg <- enrichKEGG(gene_df$ENTREZID, organism = "hsa")

write.table(as.data.frame(go),   "Results/RODAM_clusterProfiler_GO.txt",   sep="\t", quote=FALSE)
write.table(as.data.frame(kegg), "Results/RODAM_clusterProfiler_KEGG.txt", sep="\t", quote=FALSE)

############################
# END OF SCRIPT
############################
