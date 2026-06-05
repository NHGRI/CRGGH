# Load libraries ----------------------------------------------------------

library(rstatix)
library(dplyr)
library(purrr)

# MT difference per cluster ----------------------------------------------------------

md <- seurat.object@meta.data %>%
  mutate(rpca_clusters = as.factor(rpca_clusters),
         Disease = as.factor(Disease))

res <- md %>%
  group_by(rpca_clusters) %>%
  wilcox_test(percent.mt ~ Disease) %>%   # per cluster Wilcoxon
  adjust_pvalue(method = "BH") %>%        # multiple testing correction
  add_significance("p.adj")               # adds stars

res <- res %>%
  left_join(
    md %>% group_by(rpca_clusters) %>% wilcox_effsize(percent.mt ~ Disease),
    by = "rpca_clusters"
  )

# Save table --------------------------------------------------------------

write.table(
  res, 
  file = paste0(rnaProject, "-percent.mt_wilcoxTest.txt"), 
  sep = "\t", 
  row.names = F, 
  col.names = T, 
  quote = F
)
