# Usage -------------------------------------------------------------------

# Script functions:
# 1. Generate tables for IPA import
# 2. Reformat exported tables for combined analysis and plotting
# 2b. Tables are individually exported from IPA, since exporting as "Export all" does not include BH values, and the only way to do that is to export each table separately after specifying to contain the BH values

# Load setup file ---------------------------------------------------------

source("setup.R")

# Load libraries ----------------------------------------------------------

library(dplyr)
library(circlize)
library(ggalluvial)
library(tidyr)
library(fmsb)

sourceable.functions <- sourceable.functions[grep(c("ColorPalette.R"), sourceable.functions)]
invisible(sapply(sourceable.functions, source))
color.palette <- c(color.palette, color.palette)


library(Seurat)
library(dplyr)

# Generate files for IPA import ----------------------------------------------

#Create list of differentially expressed genes between TDP and TDM per cluster
seurat.object <- readRDS(paste0(rnaProject, "-analysis.object.RDS"))
for(i in levels(seurat.object$rpca_clusters)){
  cl <- subset(seurat.object, subset = rpca_clusters == i)
  Idents(cl) <- "Disease"
  clust.pos.markers <- FindAllMarkers(cl, only.pos = FALSE)
  saveRDS(clust.pos.markers, paste0("./percluster_by_disease/", rnaProject, "-Markers-cl", as.character(i), "-Disease.rds"))
  markers.pos.clust <- clust.pos.markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::select(-c(pct.1, pct.2, p_val)) %>%
    dplyr::filter(p_val_adj <= 0.05) %>%
    dplyr::rename(q_value = p_val_adj) %>% 
    dplyr::arrange(desc(abs(avg_log2FC)), .by_group = TRUE) %>% 
    select(gene, avg_log2FC, q_value, cluster)
  
  write.table(markers.pos.clust, file = paste0("./percluster_by_disease/cluster", as.character(i),".txt"), quote = F, sep = "\t", row.names = F, col.names = T)
  
  t2dlist <- markers.pos.clust %>%
    dplyr::filter(cluster == "T2D") %>% 
    ungroup() %>% 
    select(gene, avg_log2FC, q_value)
  
  ctllist <- markers.pos.clust %>%
    dplyr::filter(cluster == "CTL") %>% 
    ungroup() %>% 
    select(gene, avg_log2FC, q_value)
  
  write.table(t2dlist, file = paste0("./percluster_by_disease/T2D/cluster", as.character(i),"-T2D.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
  write.table(ctllist, file = paste0("./percluster_by_disease/CTL/cluster", as.character(i),"-CTL.txt"), quote = F, sep = "\t", row.names = F, col.names = T)
}

# Reformat exported IPA canonical pathway tables for supplemental data and graphing -------------------------------------------------

ipa.file.list <- list.files(path = "./percluster_by_disease/T2D/", pattern = "*_CP_cdb.txt", full.names = T)

for(ipa.file in ipa.file.list){
  print(basename(ipa.file))
  ipa.list <- read.csv(ipa.file, header = T, blank.lines.skip = F, sep = "\t", fill = T, skip = 2)[, 1:4]
  ipa.list["clusterID"] <- gsub(pattern = "_CP_cdb.txt", replacement = "", basename(ipa.file))
  if(ipa.file == ipa.file.list[1]){
    cp.df.dis <- ipa.list
  } else{
    cp.df.dis <- rbind(ipa.list, cp.df.dis)    
  }
}
colnames(cp.df.dis) <- c("CanonicalPway", "neglog10BH", "Ratio", "zScore", "clusterID")
cp.df.dis <- cp.df.dis %>% 
  mutate(FDR = 10^-neglog10BH) %>% 
  filter(FDR <= 0.05) %>% 
  select(c(CanonicalPway, FDR, zScore, Ratio, clusterID))
write.table(cp.df.dis, file = paste0(dirname(ipa.file.list[1]), "/", rnaProject, "-cp-FDR_cdb.txt"), col.names = T, row.names = F, sep = "\t", quote = F)

# Reformat exported IPA disease and biofunction tables for supplemental data and graphing -------------------------------------------------

ipa.file.list <- list.files(path = "./percluster_by_disease/T2D/", pattern = "*_BF_cdb.txt", full.names = T)

for(ipa.file in ipa.file.list){
  print(basename(ipa.file))
  ipa.list <- read.csv(ipa.file, header = T, blank.lines.skip = F, sep = "\t", fill = T, skip = 2)[, c(2:8, 10)]
  ipa.list["clusterID"] <- gsub(pattern = "_BF_cdb.txt", replacement = "", basename(ipa.file))
  if(ipa.file == ipa.file.list[1]){
    bf.df.dis <- ipa.list
  } else{
    bf.df.dis <- rbind(ipa.list, bf.df.dis)    
  }
}
colnames(bf.df.dis) <- c(
  "Functions",
  "pvalue",
  "FDR",
  "PredictedActivationState",
  "zScore",
  "Notes",
  "BiasCorrectedZScore",
  "No.Molecules",
  "clusterID")

bf.df.dis <- bf.df.dis %>% 
  filter(FDR <= 0.05) %>% 
  mutate(
    zScore = case_when(
      Notes == "bias" ~ BiasCorrectedZScore,
      TRUE ~ zScore
    )
  ) %>% 
  select(c(Functions, FDR, zScore, No.Molecules, clusterID))
write.table(bf.df.dis, file = paste0(dirname(ipa.file.list[1]), "/", rnaProject, "-bf-FDR_cdb.txt"), col.names = T, row.names = F, sep = "\t", quote = F)


# Generate barplot table --------------------------------------------------

#Read data
cp.cdb <- read.csv(file = paste0("./percluster_by_disease/", rnaProject, "-cp-FDR_cdb.txt"), header = T, sep = "\t")
cp.cdb <- cp.cdb %>% 
  filter(abs(zScore) >= 2) %>%
  select(CanonicalPway, zScore, FDR, clusterID)
# bf files generated with fdr <= 0.05 threshold and bias-corrected zScores inserted where appropriate
bf.cdb <- read.csv(file = paste0("./percluster_by_disease/", rnaProject, "-bf-FDR_cdb.txt"), sep = "\t", header = T)
bf.cdb <- bf.cdb %>% 
  filter(abs(zScore) >= 2 ) %>% 
  select(Functions, zScore, FDR, clusterID)
colnames(cp.cdb) <- colnames(bf.cdb)
fdr.cdb <- rbind(cp.cdb, bf.cdb)

write.table(fdr.cdb, file = paste0("./", rnaProject, "-FDR_cdb.txt"), col.names = T, row.names = F, sep = "\t", quote = F)









