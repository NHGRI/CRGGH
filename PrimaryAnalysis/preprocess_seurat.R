# Notes -------------------------------------------------------------------

# For generating initial, unanalyzed seurat object

# Load setup file ---------------------------------------------------------

source("setup.R")

##session info
session.file <- paste0(rnaProject, "-preprocess_seurat_sessionInfo.txt")
writeLines(capture.output(sessionInfo()), session.file)

# Load data ---------------------------------------------------------------

## load data list
sc.data <- sapply(list.dirs(path = path_to_data, recursive = FALSE, full.names = TRUE), 
                  basename, 
                  USE.NAMES = TRUE)

# Select data based on inclusion criteria
metadata <- read.table(file = metadata, header = TRUE, sep = "\t", row.names = 1)


# Create seurat object -------------------------------------------------------------

object.list <- foreach(i=1:length(sc.data), .combine="c", .packages = c('Seurat', 'DoubletFinder', 'ggplot2', 'dplyr')) %dopar% {
  object.item <- Read10X_h5(paste0(names(sc.data)[i], "/outs/cb-seurat_feature_bc_matrix_filtered.h5"))
  object.item <- CreateSeuratObject(object.item, 
                                    project = sc.data[i], 
                                    min.cells = min.cells, 
                                    min.features = min.features)
  object.item$orig.ident <- as.factor(sc.data[[i]])
  object.item <- AssignMetadata(metadata.df = metadata, seurat.object = object.item)
  object.item <- PercentageFeatureSet(object.item, pattern = "MT-", col.name = "percent.mt")
  object.item <- subset(object.item, 
                        subset = nFeature_RNA >= min.features &
                          nFeature_RNA <= 2500 & 
                          percent.mt <= 10)
  if(length(Cells(object.item)) < 100){
    statement <- print(paste0("Excluding ", sc.data[[i]], " because cell count = ", length(Cells(object.item))))
    print(statement)
    write(statement, file = session.file, append = TRUE)
    return(NULL)
  } else{
    print(paste("analyzing", sc.data[[i]]))
    DefaultAssay(object.item) <- "RNA"
    object.item <- NormalizeData(object.item, assay = "RNA")
    object.item <- FindVariableFeatures(object.item, selection.method = "vst", nfeatures = 3000, assay = "RNA")
    object.item <- ScaleData(object.item, assay = "RNA")
    object.item <- RunPCA(object.item, assay = "RNA")
    tot.var <- percent.variance(object.item@reductions$pca@stdev, plot.var = FALSE, return.val = TRUE)
    cluster.dims <- 0
    if(cum.var.thresh > 0){
    	cluster.dims <- length(which(cumsum(tot.var) <= cum.var.thresh))
    }
    dim.statement <- paste("Using", cluster.dims, "dims to account for", cum.var.thresh, "variance in", names(sc.data)[i])
    print(dim.statement)
    write(dim.statement,file=paste0(rnaProject, "_sessionInfo.txt"),append=TRUE)
    
    object.item <- runDoubletFinder(object.item, sctransformed = FALSE, predicted.doubletRate = predicted.doubletRate, cluster.dims = cluster.dims)
    return(object.item)
  }
}

names(object.list) <- sapply(object.list, FUN = function(x){levels(x$orig.ident)})
saveRDS(object.list, file = paste0(rnaProject, "-raw_object.list.RDS"))

write("Saved raw object list", file = session.file, append = TRUE)



# Renormalize without duplicates ------------------------------------------

object.list <- foreach(i=1:length(object.list), .combine="c", .packages = c('Seurat', 'DoubletFinder', 'ggplot2', 'dplyr')) %dopar% {
		object.item <- object.list[[i]]
	  object.item <- subset(object.item, 
                        subset = nFeature_RNA >= min.features &
                          nFeature_RNA <= 2500 & 
                          percent.mt <= 10)
  if(length(Cells(object.item)) < 100){
    statement <- print(paste0("Excluding ", sc.data[[i]], " because cell count = ", length(Cells(object.item))))
    print(statement)
    write(statement, file = session.file, append = TRUE)
    return(NULL)
  } else{
    print(paste("adding", names(object.list)[i]))
    return(object.item)
  }
}

names(object.list) <- sapply(object.list, FUN = function(x){levels(x$orig.ident)})
saveRDS(object.list, file = paste0(rnaProject, "object.list.RDS"))

write("Saved object list without duplicates", file = session.file, append = TRUE)


# Merge object ------------------------------------------------------------


object.list <- lapply(1:length(object.list), FUN = function(x){
	object.list[[x]]@assays$RNA@layers$scale.data <- NULL
	return(object.list[[x]])
})

seurat.object <- merge(object.list[[1]], y = object.list[2:length(object.list)], add.cell.ids = names(object.list))
write("Saved merged seurat object", file = session.file, append = TRUE)

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
saveRDS(seurat.object, file = paste0(rnaProject, "-FullData-integrated.RDS"))


plot.title <- "Cluster0_highlight"
cells.highlight <- WhichCells(seurat.old, idents = 0)
p <- DimPlot(
  seurat.object,
  reduction = plot.reduction, 
  cols = rep("black", length(unique(seurat.old@active.ident))),
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
png(filename = paste0(rnaProject, "-FullData-", plot.title, ".png"), height = 1000, width = 1200)
plot(p)
dev.off()




