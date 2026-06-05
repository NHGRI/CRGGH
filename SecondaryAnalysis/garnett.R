# Load libraries ----------------------------------------------------------

library(garnett)
library(org.Hs.eg.db)

# Muscle ------------------------------------------------------------------

# load classifier
muscle.classifier <- readRDS("muscle_classifier.RDS")

# load CDS
cds <- readRDS(paste0(rnaProject, "-CDS_learned-graph.RDS"))

#transfer umap graph
graph <- principal_graph(cds)$UMAP
class(graph)
upgraded_graph <- igraph::upgrade_graph(graph)
principal_graph(cds)$UMAP <- upgraded_graph

# classify cells
cds <- classify_cells(
  cds = cds, 
  classifier = muscle.classifier,
  db = org.Hs.eg.db,
  cluster_extend = TRUE,
  cds_gene_id_type = "SYMBOL")
saveRDS(cds, file = paste0(rnaProject, "-cds-GARNETT.RDS"))
head(pData(cds))
garnett.df <- as.data.frame(colData(cds))

# Save results ------------------------------------------------------------

#load seurat object
seurat.object <- readRDS(paste0(rnaProject, "-integrated.RDS"))

#add metadata and save results
seurat.object <- AddMetaData(seurat.object, metadata = garnett.df, col.name = c("cluster_ext_type", "cell_type"))
saveRDS(seurat.object, paste0(rnaProject, "-analysis.object.RDS"))


