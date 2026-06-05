# Load libraries ----------------------------------------------------------

library(enrichR)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(ComplexHeatmap)
library(viridis)


source("enrichR_functions.R")

# Generate input tables for enrichR -----------------------------------------------------------------

all.markers <- readRDS(paste0(rnaProject, "-allMarkers.rds"))
pos.markers <- readRDS(paste0(rnaProject, "-posMarkers.rds"))

coi <- c(2, 4, 6, 8, 9, 12, 14, 15, 18, 20)

all.enr <- all.markers %>% 
  filter(cluster %in% coi, abs(avg_log2FC) >= 1.5, p_val_adj <= 0.05) %>% 
  select(avg_log2FC, cluster, gene) %>% 
  group_by(cluster) %>% 
  top_frac(., .1, abs(avg_log2FC))

for(i in all.enr$cluster){
  df <- all.enr %>% filter(cluster == i)
  write.table(
    df$gene, 
    file = paste0("EnrichR/Cluster", as.character(i), "_allEnrichR-10pct.txt"), 
    row.names = F, 
    col.names = F,
    quote = F)
}


pos.enr <- pos.markers %>% 
  filter(cluster %in% coi, abs(avg_log2FC) >= 1.5, p_val_adj <= 0.05) %>% 
  select(avg_log2FC, cluster, gene) %>% 
  group_by(cluster) %>% 
  top_frac(., .1, abs(avg_log2FC))
for(i in pos.enr$cluster){
  df <- pos.enr %>% filter(cluster == i)
  write.table(
    df$gene, 
    file = paste0("EnrichR/Cluster", as.character(i), "_posEnrichR-10pct.txt"),
    row.names = F, 
    col.names = F,
    quote = F)
}

# statup enrichR ----------------------------------------------------------

websiteLive <- getOption("enrichR.live")
# if (websiteLive) {
#   listEnrichrSites()
#   setEnrichrSite("Enrichr") # Human genes
# }
# all.markers <- readRDS(paste0(rnaProject, "-allMarkers30dims_noMAST.rds"))

db.list <- "GWAS_Catalog_2025"
coi <- c(2, 4, 6, 8, 9, 12, 14, 15, 18, 20)
celltype.list <- list(
  "vascular" = c(6, 9, 12),
  "neuronal" = c(8, 15),
  "myocyte_progenitor" = c(2, 4, 18),
  "immune" = c(14, 20)
)

gene.prop <- 0.05
enrichr_adjPval = 1

# EnrichR analysis with GWAS_Catalog_2025 -------------------------------------------------------------


db.list <- "GWAS_Catalog_2025"

coi <- c(2, 4, 6, 8, 9, 12, 14, 15, 18, 20)
coi <- paste0("g", coi)
celltype.list <- list(
  "vascular" = c("g6", "g9", "g12"),
  "neuronal" = c("g8", "g15"),
  "myocyte_progenitor" = c("g2", "g4", "g18"),
  "immune" = c("g14", "g20")
)

bulk.xprsn <- read.table(
  file = paste0(rnaProject, "-pseudoBulk-AvgExpression.txt"), 
  sep = "\t", 
  header = T, 
  row.names = 1)


enrichr_adjPval = 1
count_n <- 150

file_suffix <- "pseudo_top150"
bulk.long <- bulk.xprsn %>% 
  select(any_of(coi)) %>% 
  tibble::rownames_to_column(., "gene") %>% 
  pivot_longer(!gene, names_to = "cluster", values_to = "pseudobulk") %>% 
  group_by(cluster) %>% 
  slice_max(order_by = pseudobulk, n = count_n)

#Generate and save complete list
bulk.enr <- enrichr_run(
  markers = bulk.long, 
  cluster.list = coi, 
  db.list = db.list, 
  gene.prop = gene.prop, 
  enrichr_adjPval = enrichr_adjPval,
  fc.threshold = fc.threshold, 
  rnaProject = rnaProject, 
  enrichr_dir = "EnrichR/enrichr_pseudobulk_gwas2025/", 
  file_suffix = file_suffix)

gwas_2025 <- c(
  "Coronary Artery Calcified Atherosclerotic Plaque (90 Or 130 Hu Threshold) In T2D",
  "Coronary Artery Calcified Atherosclerotic Plaque (130 Hu Threshold) In T2D",
  "Diabetic Retinopathy",
  "Fasting Blood Insulin (Bmi Interaction)",
  "Fasting Glucose",
  "Glucose Homeostasis Traits",
  "Glycemic Traits (Pleiotropy)",
  "Homeostasis Model Assessment Of Beta-Cell Function (Dietary Factor Interaction)",
  "Homeostasis Model Assessment Of Insulin Resistance",
  "Metabolite Levels (Dihydroxy Docosatrienoic Acid)",
  "Metabolite Levels (Pyroglutamine)",
  "Metabolite Risk Score For Predicting Weight Gain",
  "Mild Age-Related T2D",
  "T2D Nephropathy",
  "Total Cholesterol Change In Response To Fenofibrate In Statin-Treated T2D",
  "Type 2 Diabetes")


pb_gwas_2025 <- pb_table %>% 
  filter(Term %in% gwas_2025, db == "GWAS_Catalog_2025") %>% 
  select(Term, Combined.Score, cluster) %>% 
  pivot_wider(names_from = cluster, values_from = Combined.Score) %>% 
  mutate(across(where(is.numeric), ~sprintf("%.2f", .x)))

write.table(pb_table, file = "EnrichR/enrichr_pseudobulk_gwas2025-selectedTerms.txt", sep = "\t", col.names = T, row.names = F)


