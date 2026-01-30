############################################
# Project: Adiponectin EWAS (RODAM & AADM)
# Analysis: Pathway enrichment (missMethyl)
# Cohort: RODAM
# Author: Muhulo Muhau Mungamba
############################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

library(missMethyl)
library(dplyr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# ---- Load RODAM EWAS results ----
rodam_toptable <- read.delim(
  "TopTables/2025.10.14_TopTable_RODAM_Adiponectin_combined_T2Dcasescontrols_n315_age_sex_cells_pos_batcharry_BMI_site_DM.txt",
  sep = "\t",
  header = TRUE
)

rodam_sorted <- rodam_toptable %>%
  filter(!is.na(CpG)) %>%
  arrange(P.Value.Mval)

all_cpgs <- rodam_sorted$CpG
top1000  <- rodam_sorted$CpG[1:1000]
top5000  <- rodam_sorted$CpG[1:5000]

# ---- GO enrichment ----
go_1000 <- gometh(top1000, all_cpgs, "GO", array.type = "450K", prior.prob = TRUE)
go_5000 <- gometh(top5000, all_cpgs, "GO", array.type = "450K", prior.prob = TRUE)

# ---- KEGG enrichment ----
kegg_1000 <- gometh(top1000, all_cpgs, "KEGG", array.type = "450K", prior.prob = TRUE)
kegg_5000 <- gometh(top5000, all_cpgs, "KEGG", array.type = "450K", prior.prob = TRUE)

# ---- Save ----
dir.create("Results/RODAM_missMethyl", recursive = TRUE, showWarnings = FALSE)

write.table(go_1000, "Results/RODAM_missMethyl/GO_top1000.txt", sep = "\t", quote = FALSE)
write.table(go_5000, "Results/RODAM_missMethyl/GO_top5000.txt", sep = "\t", quote = FALSE)
write.table(kegg_1000, "Results/RODAM_missMethyl/KEGG_top1000.txt", sep = "\t", quote = FALSE)
write.table(kegg_5000, "Results/RODAM_missMethyl/KEGG_top5000.txt", sep = "\t", quote = FALSE)

############################
# END OF SCRIPT
############################
