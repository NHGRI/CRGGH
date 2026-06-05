# Load libraries ----------------------------------------------------------

library(dplyr)
library(Seurat)
library(ggplot2)

# create output files
try(setwd(rna.dir), silent = TRUE)

# Load object -------------------------------------------------------------

seurat.object <- readRDS(paste0(rnaProject , "-analysis.object.RDS"))


# Pseudobulk expression across clusters ---------------------------------------------------

avg.xprsn <- AggregateExpression(seurat.object)
avg.xprsn <- avg.xprsn$RNA
avg.xprsn <- as.data.frame(avg.xprsn)
avg.xprsn <- avg.xprsn %>% filter_all(any_vars(is.numeric(.) & . > 0))
write.table(avg.xprsn, file = paste0(rnaProject, "-pseudoBulk-AvgExpression.txt"), 
            sep = "\t", 
            quote = F, 
            col.names = T, 
            row.names = T
)  

# Pseudobulk differential expression by disease cluster ------------------------------------------------

pseudo <- AggregateExpression(
  object = seurat.object, 
  assays = "RNA",
  return.seurat = T,
  group.by = c("rpca_clusters", "Disease", "orig.ident")
)
pseudo$deseq2vars <- paste(pseudo$rpca_clusters, pseudo$Disease, sep = "_")

# DESeq2 ------------------------------------------------------------------

Idents(pseudo) <- "deseq2vars"

markers.table <- openxlsx::createWorkbook()
for(clust in unique(pseudo$rpca_clusters)){

  id1 <- paste0(clust, "_T2D")
  id2 <- paste0(clust, "_CTL")
  print(paste("FindMarkers", id1, "vs", id2))

  pseudo.deseq <- FindMarkers(
    object = pseudo,
    ident.1 = id1,
    ident.2 = id2,
    test.use = "DESeq2"
  )
  markers.deseq <- pseudo.deseq %>%
    select(-c(pct.1, pct.2, p_val)) %>%
    arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
  
  saveRDS(
    markers.deseq, 
    file = paste0("pseudobulk/", rnaProject, "-pseudodeseq_", clust, "_T2DvsCTL.RDS")
    )
  print(min(markers.deseq$p_val_adj, na.rm = T))
   
  openxlsx::addWorksheet(markers.table, sheetName = clust)
  openxlsx::writeData(markers.table, sheet = clust, x = markers.deseq, startCol = 1, startRow = 1, colNames = TRUE, rowNames = TRUE)
  
  print(paste("Finished saving", id1, "vs", id2))
}

##save workbook
openxlsx::saveWorkbook(wb = markers.table, file = paste0("pseudobulk/", rnaProject, "_pseudodeseq_T2DvsCTL.xlsx"), returnValue = TRUE)



