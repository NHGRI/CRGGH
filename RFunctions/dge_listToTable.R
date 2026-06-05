# Load libraries ----------------------------------------------------------

library(dplyr)

# Load DGE lists ----------------------------------------------------------

cluster.dir <- '.'
cluster.files <- list.files(path = cluster.dir, pattern = "cluster.*.txt")

for(file.id in cluster.files){
  cluster.id <- gsub(pattern = "cluster(.*?)-T2D.txt", "\\1", file.id)
  df <- read.table(
    file = file.path(cluster.dir, file.id), 
    header = T, 
    sep = "\t"
  )
  df$cluster <- cluster.id
  if(cluster.files[1] == file.id){
    cluster.df <- df
  } else {
    cluster.df <- rbind(cluster.df, df)
  }
}

# Filter list -------------------------------------------------------------

cluster.df <- cluster.df %>% 
  filter(q_value <= 0.05) %>% 
  select(avg_log2FC, q_value, cluster, gene)

# Write out single table --------------------------------------------------

write.table(
  x = cluster.df, 
  file = "allMarkers-percluster_by_disease.txt", 
  sep = "\t", 
  row.names = F, 
  col.names = T, 
  quote = F
)


