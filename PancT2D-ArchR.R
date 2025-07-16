# Set up ------------------------------------------------------------------

atacProject <- "PancT2D-atac"
res <- 0.5
make.arrow.files <- FALSE 
excluded.samples <- NULL

# Directories -------------------------------------------------------------

working.dir <- "./"
path_to_data <- c(list.dirs("snATAC/", full.names = TRUE, recursive = FALSE))
records.dir <- "ATACAnalysis/"
path_to_arrow_files <- "ArrowFiles/"
metadata.location <- "./"
sourceable.functions <- list.files(
  path = "RFunctions/",
  pattern = "*.R$", full.names = TRUE)
nthreads <- ceiling(future::availableCores()/2)

# Libraries ---------------------------------------------------------------

library(dplyr)
library(ArchR)
library(harmony)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Seurat)
library(SeuratObject)


addArchRThreads(threads = nthreads, force = TRUE)
addArchRGenome("hg38")

# load local functions
sourceable.functions <- sourceable.functions[grep(c("ColorPalette.R"), sourceable.functions)]
invisible(sapply(sourceable.functions, source))

# create output files
try(setwd(working.dir), silent = TRUE)

##session info
writeLines(capture.output(sessionInfo()), paste0(atacProject, "_sessionInfo.txt"))



# Generate arrowFiles -----------------------------------------------------

#make arrow files for everything the first time through
if(make.arrow.files == TRUE){
	for(i in path_to_data){
		ifelse(file.exists(paste0(i, "/outs")),"", path_to_data <-path_to_data[!path_to_data %in% i])
	}
	names(path_to_data) <- sapply(path_to_data, basename)
	arrowfiles <- createArrowFiles(
		inputFiles = paste0(path_to_data, "/outs/fragments.tsv.gz"),
		sampleNames = names(path_to_data),
		minTSS = 4, 
		minFrags = 1000, 
		addTileMat = TRUE,
		addGeneScoreMat = TRUE, force = TRUE
	)
}

# Subset Arrow Files ------------------------------------------------------

arrowfiles <- list.files(path = path_to_arrow_files, pattern = "*.arrow")
arrowfiles <- sort(arrowfiles)

# Select data based on inclusion criteria
metadata <- read.table(file = paste0(metadata.location, "HPAPMetaData.txt"), header = TRUE, sep = "\t")

# clean up metadata table
metadata <- metadata[,1:12]


# Create arch.proj -----------------------------------------------------

addArchRThreads(threads = 1) # MUST set threads=1 for addDoubletScores to work
dbltScores <- addDoubletScores(input = (arrowfiles), k = 10, knnMethod = "UMAP", LSIMethod = 1)
arch.proj <- ArchRProject(ArrowFiles = sapply(arrowfiles, function(x){paste0(path_to_arrow_files, x)}), outputDirectory = working.dir, copyArrows = TRUE)

# remove samples with unreasonable numbers of cells
for(i in getSampleNames(arch.proj)){
  # print(paste(i, "-", length(which(arch.proj$Sample %in% i))))
  if(length(which(arch.proj$Sample %in% i)) < 500){
    excluded.samples <- c(excluded.samples, i)
  }
  if(length(which(arch.proj$Sample %in% i)) > 20000){
    excluded.samples <- c(excluded.samples, i)
  }
}
arch.proj <- arch.proj[which(arch.proj$Sample %ni% excluded.samples),]


saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)


# QC plots for downstream filtering ---------------------------------------

cellcoldata <- getCellColData(arch.proj, select = c("log10(nFrags)", "TSSEnrichment"))

# plot unique nuclear fragments vs TSS enrichment score
ggPoint(
	x = cellcoldata[,1],
	y = cellcoldata[,2],
	colorDensity = TRUE,
	continuousSet = "sambaNight", 
	xlabel = "Log10 Unique Fragments",
	ylabel = "TSS Enrichment",
	xlim = c(log10(500), quantile(cellcoldata[,1], probs = 0.99)),
	ylim = c(0, quantile(cellcoldata[,2], probs = 0.99))) +
	geom_hline(yintercept = 4, lty = "dashed") +
	geom_vline(xintercept = 3, lty = "dashed")


# plot TSSEnrichment ridge plot per sample
plotGroups(
	ArchRProj = arch.proj,
	groupBy = "Sample", 
	colorBy = "cellColData",
	name = "TSSEnrichment", 
	plotAs = "ridges"
) +
	geom_vline(xintercept = 8, lty = "dashed") 

# -- TSS plots show some samples have 2 TSS enrichment at 2 different scores. Will set cutoff @ 8 to maintain a single enriched peak.

# plot nFrags ridge plot per sample
plotGroups(
	ArchRProj = arch.proj,
	groupBy = "Sample", 
	colorBy = "cellColData",
	name = "log10(nFrags)", 
	plotAs = "ridges") +
	geom_vline(xintercept = log10(1000), lty = "dashed") +
	geom_vline(xintercept = log10(30000), lty = "dashed")

arch.proj@cellColData[,names(metadata)] <- lapply(names(metadata), function(x){
	arch.proj@cellColData[[x]] <- metadata[match(vapply(strsplit(as.character(arch.proj$Sample), "_"), `[`, 1, FUN.VALUE = character(1)), metadata$DonorID), x]
}
)



# Filter  ---------------------------------------------------------

arch.proj <- filterDoublets(ArchRProj = arch.proj, cutEnrich = 1, filterRatio = 1.5) # see notes
arch.proj <- arch.proj[which(arch.proj$TSSEnrichment > 8 & 
														 	arch.proj$nFrags > 1000 & # >1000 : https://www.nature.com/articles/s41586-021-03604-1#Sec9 
														 	arch.proj$nFrags<40000 &
														 	arch.proj$BlacklistRatio < 0.03)]

arch.proj <- addIterativeLSI(
	ArchRProj = arch.proj,
	useMatrix = "TileMatrix", 
	name = "IterativeLSI", 
	iterations = 10, 
	clusterParams = list( #See Seurat::FindClusters
		resolution = c(0.3), 
		sampleCells = 10000, 
		n.start = 10
	), 
	varFeatures = 25000, 
	dimsToUse = 1:30,
	force = TRUE)

# Run Harmony -------------------------------------------------------------

# factorize regression columns
arch.proj$SampleAge <- as.factor(arch.proj$SampleAge)
arch.proj <- addHarmony(ArchRProj = arch.proj, 
												reducedDims = "IterativeLSI", 
												name = "Harmony", 
												groupBy = c("Sample", "SampleSex", "SampleAge", "DonorID"), 
												max.iter.harmony = 20, #did not converge after 10
												force = TRUE) # addHarmony "groupby" defines variables to correct for

# Add UMAP and cluster ----------------------------------------------------------------

arch.proj <- addUMAP(ArchRProj = arch.proj,
										 reducedDims = "Harmony",
										 name = "UMAP_harmony",
										 nNeighbors = 30,
										 minDist = 0.5,
										 metric = "cosine",
										 force = TRUE)
arch.proj <- addClusters(input = arch.proj, reducedDims = "Harmony", method = "Seurat", name = paste0("Harmony_res", as.character(res)), resolution = res, force = TRUE)


# UMAP and heatmap plots --------------------------------------------------

plotEmbedding(ArchRProj = arch.proj, colorBy = "cellColData", name = paste0("Harmony_res", as.character(res)), embedding = "UMAP_harmony", ) +
  theme_ArchR(legendTextSize = 12)

table(getCellColData(ArchRProj = arch.proj, select = paste0("Harmony_res", as.character(res))))
cM <- confusionMatrix(paste0(arch.proj$Harmony_res0.5), paste0(arch.proj$SampleEthnicity)) # Could not automate this line
cM <- cM / Matrix::rowSums(cM)
pheatmap::pheatmap(
	mat = as.matrix(cM),
	color = paletteContinuous("whiteBlue"),
	border_color = "black"
)

cM <- confusionMatrix(paste0(arch.proj@cellColData[,paste0("Harmony_res", as.character(res))]), paste0(arch.proj$SampleEthnicity))
cM <- cM / Matrix::rowSums(cM)
pheatmap::pheatmap(
	mat = as.matrix(cM),
	color = paletteContinuous("whiteBlue"),
	border_color = "black"
)

# Identify marker genes ---------------------------------------------------

markergenes <- getMarkerFeatures(arch.proj, groupBy = paste0("Harmony_res", as.character(res)), useMatrix = "GeneScoreMatrix", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
markerList <- getMarkers(markergenes, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
saveRDS(markergenes, file = paste0(atacProject, "-markergenes.RDS"))

# Calling peaks -----------------------------------------------------------

pathToMacs2 <- findMacs2()
BSgenome.Hsapiens.UCSC.hg18 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

arch.proj <- addGroupCoverages(arch.proj, groupBy = paste0("Harmony_res", as.character(res)))
arch.proj <- addReproduciblePeakSet(arch.proj, groupBy = paste0("Harmony_res", as.character(res)), pathToMacs2 = pathToMacs2)
arch.proj <- addPeakMatrix(arch.proj)


# Peak annotation ---------------------------------------------------------

markerPeaks <- getMarkerFeatures(arch.proj, groupBy = paste0("Harmony_res", as.character(res)), useMatrix = "PeakMatrix", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
saveRDS(markerPeaks, file = paste0(atacProject, "-MarkerPeaks.RDS"))

markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
heatmapPeaks <- plotMarkerHeatmap(seMarker = markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5", transpose = TRUE)
plot(heatmapPeaks)
arch.proj <- addMotifAnnotations(arch.proj, motifSet = "encode", annoName = "encode", force = TRUE)
arch.proj <- addArchRAnnotations(ArchRProj = arch.proj, collection = "EncodeTFBS")
enrichEncode <- peakAnnoEnrichment(seMarker = markerPeaks, ArchRProj = arch.proj, peakAnnotation = "EncodeTFBS", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")


#Plotting Encode TFBS
heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 5, transpose = FALSE, returnMatrix = FALSE)
tfbshm <-ComplexHeatmap::draw(heatmapEncode, heatmap_legend_side = "bot", annotation_legend_side = "bot")

png(filename = "ATAC_tfbsHM-raw.png", height = 800, width = 800)
plot(tfbshm)
dev.off()

heatmapEncode.matrix <- plotEnrichHeatmap(enrichEncode, n = 20, transpose = FALSE, returnMatrix = TRUE)
write.table(heatmapEncode.matrix, file = paste0(atacProject, "-heatmapEncode_matrix.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

# Motif Deviations --------------------------------------------------------

arch.proj <- addBgdPeaks(arch.proj)
if("cisbMotif" %ni% names(arch.proj@peakAnnotation)){
	arch.proj <- addMotifAnnotations(arch.proj, motifSet = "cisbp", annoName = "cisbpMotif")
}
if("encodeMotif" %ni% names(arch.proj@peakAnnotation)){
	arch.proj <- addMotifAnnotations(arch.proj, motifSet = "encode", annoName = "encodeMotif")
}
arch.proj <- addDeviationsMatrix(arch.proj, peakAnnotation = "encodeMotif", force = TRUE)

# Deviant Motifs ----------------------------------------------------------

seGroupMotif <- getGroupSE(arch.proj, useMatrix = "encodeMotifMatrix", groupBy = paste0("Harmony_res", as.character(res)))
saveRDS(seGroupMotif, file = paste0(atacProject, "_encodeMotifMatrix.RDS"))

corGSM_MM <- correlateMatrices(arch.proj, useMatrix1 = "GeneScoreMatrix", useMatrix2 = "encodeMotifMatrix", reducedDims = "Harmony")
saveRDS(corGSM_MM, file = paste0(atacProject, "_corGSM_encodeMotifMatrix.RDS"))

# Footprinting ------------------------------------------------------------

motifPositions <- getPositions(arch.proj)
saveRDS(motifPositions, file = paste0(working.dir, atacProject, "-encodeMotifPositions.rds"))

# Integrate scRNA object (Seurat) -----------------------------------------

seurat.object <- readRDS("PancT2D-18dims.RDS")
DefaultAssay(seurat.object) <- "RNA"

#Assign cluster.ids to seurat.object cells
cluster.ids <- c(
  "Alpha1", #SCT0
  "Exocrine1", #SCT1
  "Beta1", #SCT2
  "Exocrine2", #SCT3--REG1A
  "Epithelial", #SCT4
  "Mixed1", #SCT5
  "Endothelial", #SCT6--KRT18
  "Mixed2", #SCT7
  "Mesenchyme1", #SCT8
  "Alpha2", #SCT9
  "Mixed3", #SCT10
  "Mesenchyme2", #SCT11
  "MastCells", #SCT12
  "Immune1", #SCT13
  "Macrophages", #SCT14
  "Immune2", #SCT15
  "Alpha3" #SCT16
)

names(cluster.ids) <- levels(seurat.object)
seurat.object <- RenameIdents(seurat.object, cluster.ids)
seurat.object$cell.ids <- seurat.object@active.ident

my_seRNA <- seurat.object
my_seRNA <- JoinLayers(my_seRNA)
my_seRNA <- as.SingleCellExperiment(my_seRNA, assay = "RNA") # ArchR runs find variable genes, which is not supported on integrated assays

arch.proj <- mod_addGeneIntegrationMatrix(  # step takes ~95min
  ArchRProj = arch.proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
  seRNA = my_seRNA,
  addToArrow = FALSE,
  groupRNA = "cell.ids",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

atac.rna.table <- table(arch.proj$Harmony_res0.5,
                        arch.proj$predictedGroup_Un
)
write.table(atac.rna.table, file = paste0(atacProject, "_getCellColData_integratedRNA.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

# constrained integration

cM <- as.matrix(confusionMatrix(arch.proj$Harmony_res0.5, arch.proj$predictedGroup_Un))
remapClust <- colnames(cM)[apply(cM, 1 , which.max)]
names(remapClust) <- rownames(cM)
arch.proj$pop.id <- mapLabels(arch.proj$Harmony_res0.5, newLabels = remapClust, oldLabels = names(remapClust))



# Save project ------------------------------------------------------------

saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)
