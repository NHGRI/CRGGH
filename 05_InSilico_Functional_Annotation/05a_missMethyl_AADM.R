############################################
# Project: Adiponectin EWAS (RODAM & AADM)
# Analysis: Pathway enrichment (missMethyl)
# Cohort: AADM
# Author: Muhulo Muhau Mungamba
#
# Description:
# Gene Ontology (GO) and KEGG pathway enrichment
# using missMethyl, accounting for probe-number bias.
############################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

# Load libraries
library(missMethyl)
library(dplyr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# ---- Load AADM EWAS results ----
aadm_toptable <- read.delim(
  "TopTables/2024.10.22_TopTable_AADM_Adiponectin_n593_age_sex_cells_bmi_t2d.txt",
  sep = "\t",
  header = TRUE
)

# ---- Rank CpGs by p-value ----
aadm_sorted <- aadm_toptable %>%
  filter(!is.na(CpG)) %>%
  arrange(P.Value)

all_cpgs  <- aadm_sorted$CpG
top1000   <- aadm_sorted$CpG[1:1000]
top5000   <- aadm_sorted$CpG[1:5000]

# ---- GO enrichment ----
go_1000 <- gometh(
  sig.cpg = top1000,
  all.cpg = all_cpgs,
  collection = "GO",
  array.type = "EPIC",
  prior.prob = TRUE
)

go_5000 <- gometh(
  sig.cpg = top5000,
  all.cpg = all_cpgs,
  collection = "GO",
  array.type = "EPIC",
  prior.prob = TRUE
)

# ---- KEGG enrichment ----
kegg_1000 <- gometh(
  sig.cpg = top1000,
  all.cpg = all_cpgs,
  collection = "KEGG",
  array.type = "EPIC",
  prior.prob = TRUE
)

kegg_5000 <- gometh(
  sig.cpg = top5000,
  all.cpg = all_cpgs,
  collection = "KEGG",
  array.type = "EPIC",
  prior.prob = TRUE
)

# ---- Save results ----
dir.create("Results/AADM_missMethyl", recursive = TRUE, showWarnings = FALSE)

write.table(go_1000, "Results/AADM_missMethyl/GO_top1000.txt", sep = "\t", quote = FALSE)
write.table(go_5000, "Results/AADM_missMethyl/GO_top5000.txt", sep = "\t", quote = FALSE)
write.table(kegg_1000, "Results/AADM_missMethyl/KEGG_top1000.txt", sep = "\t", quote = FALSE)
write.table(kegg_5000, "Results/AADM_missMethyl/KEGG_top5000.txt", sep = "\t", quote = FALSE)

############################
# END OF SCRIPT
############################
