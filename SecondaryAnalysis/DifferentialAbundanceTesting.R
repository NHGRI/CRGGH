# scProportionTest --------------------------------------------------------

library("scProportionTest")
library(Seurat)
# Load libraries ----------------------------------------------------------


seurat.object <- readRDS(paste0(rnaProject, "-analysis.object.RDS"))


prop_test <- sc_utils(seurat.object)

scProp <- permutation_test(
  prop_test, cluster_identity = "rpca_clusters",
  sample_1 = "CTL", sample_2 = "T2D",
  sample_identity = "Disease", 
  n_permutations = 100000
)

permutation_plot(scProp, order_clusters = FALSE)

png(filename = paste0(rnaProject, "-Disease-scPropTest-100kpermutations.png"), height = 1600, width = 800)
permutation_plot(scProp, order_clusters = FALSE)
dev.off()

scProp.results <- scProp@results$permutation

write.table(x = scProp.results, file = paste0(rnaProject, "-Disease-scPropTest-100Kpermutations.txt"), quote = FALSE, row.names = F, sep = "\t", append = F)

