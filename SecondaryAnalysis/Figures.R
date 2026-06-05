# Load libraries ----------------------------------------------------------

library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyr)
library(rlang)
library(ggVennDiagram)



# load local functions
sourceable.functions <- sourceable.functions[
  grep(
    "ColorPalette.R|ipa.networks.R", 
    sourceable.functions)]
invisible(sapply(sourceable.functions, source))
color.palette <- c(color.palette, color.palette)


# create output files
try(setwd(rna.dir), silent = TRUE)

# Load object -------------------------------------------------------------

seurat.object <- readRDS(paste0(rnaProject , "-analysis.object.RDS"))

# Figure 1 ----------------------------------------------------------------
figmarkers <- c(
  'MYOD1', 'PAX7', #myoblast and/or progenitor
  'APOC1', 'APOE', 'MYF5', 'MYOG', 'MET', #satellite
  'MYH7', #type I fiber
  'MYH1', #Type IIX fiber
  'MYH2', #type IIA fiber
  'RGS5', 'ACTA2', 'NOTCH3',#pericyte
  'LUM', 'DCN', 'FBN1', 'PCOLCE2', #FAP
  'APOD', #adipocyte
  'COL3A1', 'COL1A1', 'PDGFRA', #fibroblast
  'PECAM1', 'CDH5', 'ESAM', 'FABP4', 'VWF',#endothelial
  'MYH11', #smooth muscle
  'BDNF', 'PTPRT', 'KCNIP4', 'CFTR', 'NSG2',
  'PTPRC', 'MS4A2', 'ITGAM'
)

figmarkers <- figmarkers[figmarkers %in% Features(seurat.object)]

#### Figure 1C
png(filename = paste0(rnaProject, "-DotPlot-assayDevmarkers.png"), height = 700, width = 1600)
DotPlot(seurat.object, features = figmarkers, col.min = 0, dot.min = 0.01, dot.scale = 9) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.grid = element_line(color = "gray", size = 0.2, linetype = 2))
dev.off()

#### Figure 1D
# prep data
heatmap.subset <- subset(seurat.object, downsample = 100)
heatmap.data <- FetchData(object = heatmap.subset, layer = "data", vars = figmarkers) %>% 
  filter(rowSums(across(where(is.numeric)))!=0)
heatmap.data <- scale(heatmap.data) %>% 
  as.matrix()
heatmap.subset <- subset(heatmap.subset, cells = rownames(heatmap.data))
cluster_anno<- heatmap.subset@meta.data$rpca_clusters

# set bar palette for heatmap to match annotations on UMAP
color.vector <- as.numeric(names(column_order(ht))) +1
ht.palette <- color.palette[color.vector]

# define color function
col_fun = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")) # set color scale based on 10% and 95% quantiles

#plot heatmap
png(filename = paste0(rnaProject, "-heatmap-identMarkers.png"), height = 800, width = 1200, bg = "transparent")
ht.png <- draw(Heatmap(t(heatmap.data), name = "Expression",  
                       column_split = factor(cluster_anno),
                       cluster_rows = TRUE,
                       show_row_dend = TRUE,
                       cluster_row_slices = TRUE,
                       row_title_gp = gpar(fontsize = 10),
                       row_gap = unit(1, "mm"),
                       cluster_columns = TRUE,
                       show_column_dend = F,
                       col = col_fun,
                       column_names_gp = gpar(fontsize = 10),
                       row_title_rot = 90, 
                       top_annotation = columnAnnotation(foo = anno_block(gp = gpar(fill = ht.palette))),
                       show_row_names = TRUE, show_column_names = FALSE, 
                       border_gp = gpar(col = "black", lty = 1), use_raster = TRUE))
dev.off()


# Figure 2 ------------------------------------------------------

#### Figure 2C
scProp.results <- read.table(file = paste0(rnaProject, "-Disease-scPropTest-100Kpermutations.txt"), sep = "\t", header = T)
cluster.order <- c(
  0, 19, 5,  10, 
  1, 2, 3, 4, 18, 
  8, 7, 15, 
  11, 13, 17,
  6, 12, 9, 
  14, 16, 20
)

sc.plot <- scProp.results %>% 
  select(clusters, obs_log2FD, FDR, boot_CI_2.5, boot_CI_97.5) %>% 
  mutate(clusters = factor(clusters, levels = cluster.order)) %>% 
  mutate(significance = 
           ifelse(FDR <= 0.05 & (obs_log2FD <= -0.46 | obs_log2FD >= 0.46), no = "n.s.", yes = "sig"),
  ) %>% 
  mutate(significance = 
           ifelse(FDR <= 0.05 & (obs_log2FD <= -0.58 | obs_log2FD >= 0.58), no = significance, yes = "vsig"),
  ) 


p <- ggplot(sc.plot, aes(x = clusters, y = obs_log2FD)) +
  geom_rect(mapping = aes(ymin=-0.45, ymax=-0.58, xmin=0, xmax=Inf), fill = "lightgray", alpha = 0.5) +
  geom_rect(mapping = aes(ymin=0.45, ymax=0.58, xmin=0, xmax=Inf), fill = "lightgray", alpha = 0.5) +
  geom_pointrange(aes(ymin = boot_CI_2.5, ymax = boot_CI_97.5, color = significance), size = 1.5, linewidth = 3) +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("gray", "black", "red")) +
  theme_bw() +
  coord_flip()

plot(p)
png(filename = paste0(rnaProject, "-scProportionTest_plot.png"), height = 800, width = 800)
plot(p)
dev.off()

#### Figure 1B
png(filename = paste0(rnaProject, "-umap-rpca_clusters.png"), height = 1200, width = 1200)
DimPlot(
  seurat.object, 
  reduction = "umap.rpca", 
  cols = color.palette, 
  shuffle = T, 
  label = F, 
  label.size = 10, 
  repel = T, 
  group.by = "rpca_clusters", 
  pt.size = 1, 
  raster = FALSE) + 
  NoLegend() +
  NoAxes()
dev.off()


# Figure 3 ----------------------------------------------------------------

dm.plots <- c("ctl", "t2d")
dm.subset <- c()
for(i in dm.plots){
  dm.subset <- c(dm.subset, readRDS(file = paste0(rnaProject, "-", i, ".cds-", i, "_learned-graph.RDS")))
}
names(dm.subset) <- dm.plots

traj.segments <- c()
traj.points <- c()
traj.text <- c()

for(i in 1:length(names(dm.subset))){
  query.name <- names(dm.subset)[i]
  traj.color <- color.palette[i]
  edge_df <-  readRDS(file = paste0(rnaProject, "-", query.name, "-learned_graph_DF.RDS"))
  branch_point_df <- readRDS(file = paste0(rnaProject, "-", query.name, "-learned_graph_branchpts.RDS"))
  
  
  segments <- geom_segment(aes_string(x="source_prin_graph_dim_1",
                                      y="source_prin_graph_dim_2",
                                      xend="target_prin_graph_dim_1",
                                      yend="target_prin_graph_dim_2"),
                           linewidth=2,
                           color=I(traj.color),
                           linetype="solid",
                           na.rm=TRUE,
                           data=edge_df)
  traj.segments <- c(traj.segments, segments)
  
  points <- geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                       shape = 21, stroke=I(1),
                       color="white",
                       fill=traj.color,
                       size=I(5 * 1.5),
                       na.rm=TRUE, branch_point_df)
  traj.points <- c(traj.points, points)
  text <- geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                               label="branch_point_idx"),
                    size=I(5), color="white", na.rm=TRUE,
                    branch_point_df)
  traj.text <- c(traj.text, text)
}
names(traj.segments) <- names(dm.subset)
names(traj.points) <- names(dm.subset)
names(traj.text) <- names(dm.subset)

#branch points

for(i in 1:length(traj.segments)){
  p <- DimPlot(down.sample, reduction = "umap.rpca", cols = "white", label = F, group.by = "DF.classifications", alpha = 0) & 
    xlim(c(-20, 20)) &
    ylim(c(-20, 20)) &
    NoLegend() &
    theme(legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent",
                                         color = NA))
  p <- p + 
    traj.segments[i] + 
    traj.points[i] + 
    theme(panel.background = element_blank()) + 
    ggtitle(paste0(names(traj.segments)[i], "-trajWithPoints"))
  # png(filename = paste0(rnaProject, "-", names(traj.segments)[i], "-trajWithPoints.png"), height = 1200, width = 1200, bg = "transparent")
  plot(p)
  # dev.off()
  
}

for(i in 1:length(names(dm.subset))){
  if(i == 1){
    traj.segments <- c()
  }
  query.name <- names(dm.subset)[i]
  edge_df <-  readRDS(file = paste0(rnaProject, "-", query.name, "-learned_graph_DF.RDS"))
  
  segments <- geom_segment(aes_string(x="source_prin_graph_dim_1",
                                      y="source_prin_graph_dim_2",
                                      xend="target_prin_graph_dim_1",
                                      yend="target_prin_graph_dim_2"),
                           linewidth=1,
                           color=I("black"),
                           linetype="solid",
                           na.rm=TRUE,
                           data=edge_df)
  traj.segments <- c(traj.segments, segments)
}
names(traj.segments) <- names(dm.subset)

#CTL
png(filename = paste0(rnaProject, "-CTL.down-traj.png"), height = 1200, width = 1200, bg = "transparent")
DimPlot(subset(down.sample, subset = Disease == "CTL"), reduction = "umap.rpca", cols = color.palette, label = F, group.by = "rpca_clusters", alpha = .5, split.by = "Disease") + 
  NoLegend() + 
  NoAxes() + 
  traj.segments$ctl
dev.off()

#T2D
png(filename = paste0(rnaProject, "-T2D.down-traj.png"), height = 1200, width = 1200, bg = "transparent")
DimPlot(subset(down.sample, subset = Disease == "T2D"), reduction = "umap.rpca", cols = color.palette, label = F, group.by = "rpca_clusters", alpha = .5, split.by = "Disease") + 
  NoLegend() +
  NoAxes() + 
  traj.segments$t2d
dev.off()



# Figure 4 -----------------------------------------------------------

#### Figure 4A
sk.table <- all.markers %>% 
  filter(p_val_adj <= .05, avg_log2FC >= 1.2, cluster %in% c(0, 1, 5, 10, 19)) %>% 
  select(gene, cluster)

sk.venn <-  list(
  "c0" = sk.table[sk.table[["cluster"]] == 0,]$gene,
  "c1" = sk.table[sk.table[["cluster"]] == 1,]$gene,
  "c5" = sk.table[sk.table[["cluster"]] == 5,]$gene,
  "c10" = sk.table[sk.table[["cluster"]] == 10,]$gene,
  "c19" = sk.table[sk.table[["cluster"]] == 19,]$gene
)

png(filename = paste0(rnaProject, "-skVenn_colorbycount.png"), height = 1000, width = 1000, bg = "transparent")
ggVennDiagram(sk.venn, label = "none", category.names = rep("", 5)) + 
  ggplot2::scale_fill_gradient(low="white",high = "blue")
dev.off()

#### Figure 4B
qc.list <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "pANN")
for (qc in qc.list){
  p <- ggplot(
    seurat.object@meta.data, 
    aes(
      x = rpca_clusters, 
      y = !!sym(qc), 
      fill = Disease)) +
    geom_boxplot(outlier.size = .5) + 
    theme_classic() +
    labs(y = qc)
  print(p)
  ggsave(
    filename = paste0("pngs/", rnaProject, "-boxandwhiskers-", qc, ".png"),
    height = 3,
    width = 9,
    dpi = 300
  )  
}

#### Figure 4C
cluster.col.map <- color.palette[1:21]
names(cluster.col.map) <- paste(0:20)

#Read table
fdr.cdb <- read.csv(file = paste0("./percluster_by_disease/T2Donly/", rnaProject, "-cp-FDR_cdb.txt"), header = T, sep = "\t")

#Subset table
sk.clusters <- paste0("Cluster", c(0, 19, 5, 10, 1))
sk.networks <- c(
  "Regulation of the Epithelial Mesenchymal Transition by Growth Factors Pathway",
  "Vasculogenesis", 
  "Development of Neurons",
  "Translocation of SLC2A4",
  "IL-6 Signaling",
  "Function of Skeletal Muscle",
  "Apoptosis of Vascular Smooth Muscle Cells",
  "Glucose Metabolism Disorder",
  "Glyoxylate Metabolism and Glycine Degradation",
  "Synthesis of Reactive Oxygen Species",
  "Migration of Smooth Muscle Cells",
  "Abnormal Morphology of Cardiovascular System"
)
sk.table <- fdr.cdb %>% 
  filter(grepl(paste(sk.networks, collapse = "|"), Functions, ignore.case = T) & clusterID %in% sk.clusters) %>% 
  mutate (clusterID = factor(clusterID, levels = sk.clusters)) %>% 
  mutate (Functions = factor(Functions, levels = unique(sk.table$Functions))) %>% 
  arrange((clusterID))

p <- ggplot(sk.table, aes(x = zScore, y = Functions, fill = clusterID, alpha = -log10(FDR))) +
  geom_col(position = "dodge") +  # Dodge for grouping bars side by side
  scale_alpha_continuous(range = c(0.5, 1)) +
  scale_fill_manual(values = cluster.col.map) +
  scale_x_continuous(breaks = seq(-3, 3, by = 0.5))+
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(axis.text = element_blank())

png(filename = file = paste0("./", rnaProject, "-FDR_cdb.txt"), height = 1200, width = 800, bg = "transparent")
plot(p)
dev.off()



#### Figure 4D
cp.cdb <- read.csv(file = paste0("./percluster_by_disease/T2Donly/", rnaProject, "-cp-FDR_cdb.txt"), header = T, sep = "\t")

head(cp.cdb)
cp.cdb <- cp.cdb %>% 
  filter(abs(zScore) >= 2) %>%
  select(CanonicalPway, zScore, clusterID)

# bf files generated with fdr <= 0.05 threshold and bias-corrected zScores inserted where appropriate
bf.cdb <- read.csv(file = paste0("./percluster_by_disease/T2Donly/", rnaProject, "-bf-FDR_cdb.txt"), sep = "\t", header = T)

bf.cdb <- bf.cdb %>% 
  filter(abs(zScore) >= 2 ) %>% 
  select(Functions, zScore, clusterID)

colnames(cp.cdb) <- colnames(bf.cdb)
all.cdb <- rbind(cp.cdb, bf.cdb)


radar.clusters <- c(
  0, 19, 5, 10, 1, 
  2, 3, 4, 18, 
  11, 13, 17, 
  6, 12, 9, 
  14, 16, 20,
  8, 7, 15
)
cluster.order <- paste0("Cluster", radar.clusters)
missed.cdb <- c()
for(i in ipa.networks){
  inverted.zscores <- FALSE
  if(length(grep(i, all.cdb$Functions, ignore.case = T)) > 3){
    cell.function <- all.cdb[grepl(i, all.cdb$Functions, ignore.case = T),]
    # print(cell.function$Functions[1])
    missing.clusters <- setdiff(cluster.order, cell.function$clusterID)
    missing.rows <-  tibble(
      "Functions" = cell.function[1,1],
      "zScore" = 0,
      clusterID = missing.clusters
    )
    zscores <- all.cdb %>% 
      filter(cell.function[1,1] == Functions) %>% 
      # mutate(ActivationZScore = 2^(1/ActivationZScore)) %>%
      bind_rows(., missing.rows) %>%
      tidyr::pivot_wider(., names_from = clusterID, values_from = zScore)
    zscores <- zscores[-1]
    
    if (all(zscores[zscores != 0] < 0, na.rm = TRUE)) {
      radar.limits <- c(1, min(zscores, na.rm = TRUE)*1.05)
    } else if (all(zscores[zscores != 0] > 0, na.rm = TRUE)) {
      radar.limits <- c(max(zscores, na.rm = TRUE)*1.05, 0)  # All non-zero values are positive
    } else {
      radar.limits <- c(max(zscores, na.rm = TRUE) *1.05, min(zscores, na.rm = TRUE))  # Mixed values
    }
    
    zscores <- zscores %>% 
      rbind(rep(radar.limits[1], ncol(zscores)),
            rep(radar.limits[2], ncol(zscores)),
            rep(0, ncol(zscores)),
            .)  %>% 
      mutate(across(where(is.numeric), \(x)round(x, 1)))# \(x)function defines the function
    
    
    zscores <- as.data.frame(apply(zscores, 2, as.numeric))
    caxis.labels <- round(seq(radar.limits[2], radar.limits[1], length.out = 5), 1)
    
    # Plot the radar chart
    custom_radarchart(zscores[,cluster.order[cluster.order %in% colnames(zscores)]], 
                      pfcol = c(NA,rgb(1,0,0,0.3)),
                      pcol= c(1,"red"), 
                      plty = c(1, 1), 
                      plwd = 3,
                      cex = c(0, 2),
                      axislabcol = "black", 
                      caxislabels = caxis.labels,
                      axistype = 1, seg = 4,
                      title = paste0("cdb_", i))
    
    # legend("topright", legend = clusters, col = c("red", "blue", "green", "purple"), lty = 1, lwd = 2)
    png(filename = paste0("RadarCharts/", rnaProject, "-radar-cdb-", paste(i, collapse = ""), ".png"), 
        height = 1200, width = 1000)
    custom_radarchart(zscores[,cluster.order[cluster.order %in% colnames(zscores)]], 
                      pfcol = c(NA,rgb(1,0,0,0.3)),
                      pcol= c(1,"red"), 
                      plty = c(1, 1), 
                      plwd = 3,
                      cex = c(0, 2),
                      axislabcol = "black", 
                      caxislabels = caxis.labels,
                      axistype = 1, seg = 4,
                      title = paste0("cdb_", i))
    dev.off()
  } else {
    print(paste("...couldn't print anything for", i))
    missed.cdb <- c(missed.cdb, i)
  }
  
}


# Figure 5 ----------------------------------------------------------------

pb_table <- read.table(
  file = "EnrichR/enrichr_pseudobulk_gwas2025-selectedTerms.txt", 
  sep = "\t", 
  header = T
)
pb_heatmap <- pb_table %>% 
  select(Term, Combined.Score, cluster) %>% 
  mutate(cluster = sub("g", "Cluster ", cluster)) %>% 
  filter(Term %in% gwas_2025) %>% 
  pivot_wider(names_from = cluster, values_from = Combined.Score
  )

pb_heatmap <- pb_heatmap[match(gwas_2025, pb_heatmap$Term),]
rownames(pb_heatmap) <- pb_heatmap$Term
mat <- log10(as.matrix(pb_heatmap[-1]))
max_val <- max(abs(mat), na.rm = TRUE)  # symmetric range around 0
min_val <- min(abs(mat), na.rm = TRUE)  # symmetric range around 0
col_fun <- circlize::colorRamp2(
  breaks = c(-min_val, 1, max_val),
  colors = colorspace::diverge_hsv(n = 3)  # low, white (mid), high
)

p1 <- ComplexHeatmap::Heatmap(
  log10(as.matrix(pb_heatmap[-1])), 
  col = col_fun,
  na_col = "white", 
  cluster_rows = F, 
  cluster_columns = F, 
  row_labels = stringr::str_wrap(rownames(pb_heatmap), width = 45),
  row_names_side = "left",
  row_names_max_width = unit(8.5, "cm"),
  rect_gp = gpar(col = "black", lwd = 1), 
  heatmap_legend_param = list(
    title = "Combined \n Score \n (log10)"
  ))

draw(p1)

pdf(file = "GWAS2025_log10_combinedscore_heatmap.pdf", height = 7, width = 9)
draw(p1)
dev.off()



# Supplemental Figure 1 ---------------------------------------------------

#### SFig 1A



#### SFig 1B
#highlight each cluster
for(i in unique(seurat.object@active.ident)){
  plot.title <- paste("highlight cluster", as.character(i))
  cells.highlight <- WhichCells(seurat.object, idents = i)
  p <- DimPlot(
    seurat.object,
    reduction = plot.reduction, 
    cols = rep("black", length(unique(seurat.object@active.ident))),
    shuffle = T, 
    label = F, 
    label.size = 7, 
    repel = T, 
    group.by = "rpca_clusters", 
    cells.highlight = cells.highlight,
    cols.highlight = "red",
    ncol = 1, 
    pt.size = .1, 
    raster = FALSE) + 
    NoLegend() + 
    ggtitle(plot.title)
  png(filename = paste0(rnaProject, "-highlight_Cluster", as.character(i), ".png"), height = 1000, width = 1200)
  plot(p)
  dev.off()
}





