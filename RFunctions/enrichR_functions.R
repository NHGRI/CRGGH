enrichr_table <- function(markers.df, cluster.list, gene.prop, fc.threshold){
   if(fc.threshold >0){
    xprsn.type <- "pos"
  } else if (fc.threshold < 0){
    xprsn.type <- "neg"
  } else {
    xprsn.type <- ""
  }
  enr.markers <- markers.df %>% 
    filter(cluster %in% cluster.list, abs(avg_log2FC) >= fc.threshold, p_val_adj <= 0.05) %>% 
    select(avg_log2FC, cluster, gene) %>% 
    group_by(cluster) %>% 
    slice_max(avg_log2FC, prop = gene.prop)
  return(enr.markers)
}

enrichr_run <- function(markers.df, cluster.list, db.list, gene.prop, enrichr_adjPval = 0.5, fc.threshold, rnaProject, enrichr_dir = "", file_suffix = ""){

  adjPval.character <- paste0("adjP", enrichr_adjPval)
  
  for(clust in cluster.list){
    enr.genes <- markers.df %>% 
      filter(clust == cluster) %>% 
      pull(gene)
    enr.list <- enrichr(enr.genes, db.list)
    for(i in names(enr.list)){
      partial.df <- enr.list[[i]] %>% 
        filter(Adjusted.P.value <= enrichr_adjPval) %>% 
        mutate(db = i, cluster = clust)
      if(i == names(enr.list)[1]){
        clutser.df <- partial.df
      } else{
        clutser.df <- rbind(clutser.df, partial.df)
      }
    }
    if(clust == cluster.list[1]){
      enr.df <- clutser.df
    } else {
      enr.df <- rbind(enr.df, clutser.df)
    }
    print(paste("finished", as.character(clust)))
  }
  
  file.name <- paste0(rnaProject, "-EnrichR", adjPval.character, "_", file_suffix)
  print(paste("saving", file.name))
  saveRDS(
    enr.df, file = paste0(enrichr_dir, file.name, ".RDS")
    )
  write.table(
    enr.df, 
    file = paste0(enrichr_dir, file.name, ".txt"), 
    col.names = T,
    row.names = F,
    sep = "\t", 
    quote = F
  )
  return(enr.df)
}

enrichr_graph <- function(celltype.list, enr.df, plot_n = 5, file_suffix = ""){
  
   enrichr_maxpval <- paste0("adjP", round(max(enr.df$Adjusted.P.value), 2))
 
    for(pop in names(celltype.list)){
    cluster_set <- celltype.list[[pop]]
    for(db_i in unique(enr.df$db)){
      print(paste0(pop, "-", db_i))
      db_i.terms <- enr.df %>% 
        filter(db == db_i, cluster %in% cluster_set) %>% 
        slice_max(Combined.Score, n = plot_n, with_ties = TRUE) %>% 
        pull(Term)
      db_i.plot <- enr.df %>% 
        filter(db == db_i, cluster %in% cluster_set, Term %in% db_i.terms)
      if(nrow(db_i.plot) > 0){
        if(nrow(db_i.plot) <=8){
          plot.height <- 6
        } else {
          plot.height <- 10
        }
        db_i.plot <- ggplot(db_i.plot, aes(x = factor(cluster, levels = sort(unique(cluster))), 
                                           y = Term)) +
          geom_point(aes(size = -log10(Adjusted.P.value), 
                         fill = Combined.Score), 
                     alpha = 0.75, shape = 21) +
          scale_fill_viridis_c(name = "Combined.Score") +
          scale_size_continuous(name = "-log10(Adjusted.P.value)", range = c(2, 10)) +
          xlab("Cluster ID") +
          scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
          ggtitle(label = paste0(db_i, " - ", pop), subtitle = paste("(enrichR adjP<=", enrichr_maxpval, ",", file_suffix, ")")) +
          theme_bw() +
          theme(
            panel.background = element_rect(fill='transparent'),
            plot.background = element_rect(fill='transparent', color=NA),
            legend.background = element_rect(fill='transparent'),
            legend.box.background = element_rect(fill='transparent'),
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = "gray85")
          )
        print(db_i.plot)
        ggsave(
          filename = paste0("EnrichR/", rnaProject, "-plot-", pop, "_", db_i, "-EnrichR", enrichr_maxpval, "_", file_suffix, ".png"),
          plot = db_i.plot, width = plot.height, height = plot.height, dpi = 300, bg = "transparent"
        )
      } else {
        print(paste("--No entries for", pop, "pop in", db_i))
      }
    }
  }
}
