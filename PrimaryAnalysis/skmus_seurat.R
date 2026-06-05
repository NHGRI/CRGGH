# Load setup file ---------------------------------------------------------

source("setup.R")

# Load libraries and local functions ----------------------------------------------------------

library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)

sourceable.functions <- sourceable.functions[grep(c("ColorPalette.R|PercentVariance.R|runDoubletFinder.R|AssignMetadata.R"), sourceable.functions)]
invisible(sapply(sourceable.functions, source))
color.palette <- c(color.palette, color.palette)

# Subset without cluster 0 --------------------------------------------

# Based on extensive analysis before 1/21/25, have concluded that Cluster0 corresponds to background reported in methdos paper. Going to exclude this cluster, and restart the analysis
seurat.subset <- readRDS(file = paste0(rnaProject, "-FullData-integrated.RDS"))
seurat.subset <- subset(seurat.subset, subset = rpca_clusters == 0, invert = TRUE)

DimPlot(seurat.subset, reduction = "umap.rpca", cols = color.palette, shuffle = T, label = F, label.size = 7, repel = T, group.by = "rpca_clusters")

# create clean unanalyzed seurat object -----------------------------------

seurat.subset <- JoinLayers(seurat.subset)
counts <- LayerData(seurat.subset, layer = "counts")

seurat.object <- CreateSeuratObject(counts = counts, assay = "RNA", meta.data = seurat.subset@meta.data)
saveRDS(seurat.object, file = paste0(rnaProject, "-RAW.RDS"))

# Unintegrated preprocessing -------------------------------------------------------

#read object
seurat.object <- readRDS(paste0(rnaProject, "-RAW.RDS"))
seurat.object[["RNA"]] <- split(seurat.object[["RNA"]], f = seurat.object$orig.ident)
if(any(unique(seurat.object$DF.classifications) %in% "Doublet")){
	seurat.object <- subset(seurat.object, subset = DF.classifications == "Singlet")
}

#log normalize and scale
seurat.object <- NormalizeData(seurat.object, assay = "RNA")
seurat.object <- FindVariableFeatures(seurat.object, selection.method = "vst", nfeatures = 3000, assay = "RNA")
seurat.object <- ScaleData(seurat.object, assay = "RNA")
seurat.object <- RunPCA(seurat.object, assay = "RNA")

#calculate dimensionality
ElbowPlot(seurat.object)
tot.var <- percent.variance(seurat.object@reductions$pca@stdev, plot.var = FALSE, return.val = TRUE)
cluster.dims <- 0
if(cum.var.thresh > 0){
  cluster.dims <- length(which(cumsum(tot.var) <= cum.var.thresh))
}

#clustering to confirm integration requirement
seurat.object <- FindNeighbors(seurat.object, dims = 1:cluster.dims, reduction = "pca")
seurat.object <- FindClusters(seurat.object, resolution = resolution, cluster.name = "unintegrated_clusters")
seurat.object <- RunUMAP(seurat.object, dims = 1:cluster.dims, reduction.name = "umap.unintegrated")

# Integrate ---------------------------------------------------------------

seurat.object <- IntegrateLayers(
	object = seurat.object, method = RPCAIntegration,
	orig.reduction = "pca", new.reduction = "rpca",
	verbose = TRUE
)

# Clustering --------------------------------------------------------------

seurat.object <- FindNeighbors(seurat.object, dims = 1:cluster.dims, reduction = "rpca")
seurat.object <- FindClusters(seurat.object, resolution = resolution, cluster.name = "rpca_clusters")
seurat.object <- RunUMAP(seurat.object, dims = 1:cluster.dims, reduction.name = "umap.rpca", reduction = "rpca")
seurat.object$rpca_clusters <- as.factor(seurat.object$rpca_clusters)

# Markers and distrubution tables -------------------------------------------------------------

all.markers <- FindAllMarkers(seurat.object, only.pos = FALSE)
saveRDS(all.markers, file = paste0(rnaProject, "-allMarkers.rds"))
pos.markers <- FindAllMarkers(seurat.object, only.pos = TRUE)
saveRDS(pos.markers, file = paste0(rnaProject, "-posMarkers.rds"))

markers.pos <- pos.markers %>%
	group_by(cluster) %>%
	select(-c(pct.1, pct.2, p_val)) %>%
	filter(p_val_adj <= 0.05) %>%
	arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
markers.all <- all.markers %>%
	group_by(cluster) %>%
	select(-c(pct.1, pct.2, p_val)) %>%
	filter(p_val_adj <= 0.05, avg_log2FC >=1 | avg_log2FC <= -1) %>%
	arrange(desc(abs(avg_log2FC)), .by_group = TRUE)

#distribution tables
cluster.distribution.by.orig.ident <- table(seurat.object$seurat_clusters,seurat.object$orig.ident)
cluster.distribution.by.sex <- table(seurat.object$seurat_clusters,seurat.object$Sex)
cluster.distribution.by.disease <- table(seurat.object$seurat_clusters,seurat.object$Disease)

#create xlsx tables
source("xlsx.tablefy.R")

markers.table <- openxlsx::createWorkbook()
xlsx.tablefy(workbook.table = markers.pos, sheet.name = "PosMarkers", workbook.name = markers.table, style.cols = "avg_log2FC", sort.it = T)
xlsx.tablefy(workbook.table = markers.all, sheet.name = "AllMarkers", workbook.name = markers.table, style.cols = "avg_log2FC", sort.it = T)

##write cluster.distribution.by.orig.ident to table
openxlsx::addWorksheet(markers.table, sheetName = "ClusDist-OrigID")
openxlsx::writeData(markers.table, sheet = "ClusDist-OrigID", x = cluster.distribution.by.orig.ident, startCol = 1, startRow = 1, colNames = TRUE)

##write cluster.distribution.by.sex to table
openxlsx::addWorksheet(markers.table, sheetName = "ClusDist-Sex")
openxlsx::writeData(markers.table, sheet = "ClusDist-Sex", x = cluster.distribution.by.sex, startCol = 1, startRow = 1, colNames = TRUE)

##write cluster.distribution.by.disease to table
openxlsx::addWorksheet(markers.table, sheetName = "ClusDist-Disease")
openxlsx::writeData(markers.table, sheet = "ClusDist-Disease", x = cluster.distribution.by.disease, startCol = 1, startRow = 1, colNames = TRUE)

##save workbook
openxlsx::saveWorkbook(wb = markers.table, file = paste0(rnaProject, "_seuratMarkers-", as.character(cluster.dims), "dims.xlsx"), overwrite = TRUE, returnValue = TRUE)



