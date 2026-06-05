# Trajectory --------------------------------------------------------------

library(monocle3)
library(SeuratWrappers)

# Full object - Load object -------------------------------------------------------------

cds <- readRDS(file = paste0(rnaProject, "-integrated.RDS"))
cds[["UMAP"]] <- cds[["umap.rpca"]]
cds[["umap.rpca"]] <- NULL
cds <- as.cell_data_set(cds)

# Full object - Run trajectory analysis --------------------------------------------------

cds <- cluster_cells(cds)
cds <- learn_graph(cds)
saveRDS(cds, file = paste0(rnaProject, "-CDS_learned-graph.RDS"))

# Full object - Save trajectory info ----------------------------------------------------

ica_space_df <- t(cds@principal_graph_aux$UMAP$dp_mst) %>%
  as.data.frame() %>%
  dplyr::select(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
  dplyr::mutate(sample_name = rownames(.),
                sample_state = rownames(.))
dp_mst <- cds@principal_graph[["UMAP"]]
edge_df <- dp_mst %>%
  igraph::as_data_frame() %>%
  dplyr::select(source = "from", target = "to") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select(
                       source="sample_name",
                       source_prin_graph_dim_1="prin_graph_dim_1",
                       source_prin_graph_dim_2="prin_graph_dim_2"),
                   by = "source") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select(
                       target="sample_name",
                       target_prin_graph_dim_1="prin_graph_dim_1",
                       target_prin_graph_dim_2="prin_graph_dim_2"),
                   by = "target")

mst_branch_nodes <- monocle3:::branch_nodes(cds, "UMAP")
branch_point_df <- ica_space_df %>%
  dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>%
  dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))
saveRDS(edge_df, file = paste0(rnaProject, "-learned_graph_DF.RDS"))
saveRDS(branch_point_df, file = paste0(rnaProject, "-learned_graph_branchpts.RDS"))

# Split Disease - Load object & subset --------------------------------------------------------------

#load object
cds <- readRDS(file = paste0(rnaProject, "-integrated.RDS"))
cds[["UMAP"]] <- cds[["umap.rpca"]]
cds[["umap.rpca"]] <- NULL

#subset
ctl.cds <- subset(cds, subset = Disease == "CTL")
t2d.down <- subset(cds, subset = Disease == "T2D")

# Full object - CTL -------------------------------------------------------

ctl.cds <- as.cell_data_set(ctl.cds)
ctl.cds <- cluster_cells(ctl.cds)

ctl.cds <- learn_graph(ctl.cds)
saveRDS(ctl.cds, file = paste0(rnaProject, "-ctl.cds-ctl_learned-graph.RDS"))

# save graph components for plotting separately
ica_space_df <- t(ctl.cds@principal_graph_aux$UMAP$dp_mst) %>%
  as.data.frame() %>%
  dplyr::select(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
  dplyr::mutate(sample_name = rownames(.),
                sample_state = rownames(.))
dp_mst <- ctl.cds@principal_graph[["UMAP"]]
edge_df <- dp_mst %>%
  igraph::as_data_frame() %>%
  dplyr::select(source = "from", target = "to") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select(
                       source="sample_name",
                       source_prin_graph_dim_1="prin_graph_dim_1",
                       source_prin_graph_dim_2="prin_graph_dim_2"),
                   by = "source") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select(
                       target="sample_name",
                       target_prin_graph_dim_1="prin_graph_dim_1",
                       target_prin_graph_dim_2="prin_graph_dim_2"),
                   by = "target")

mst_branch_nodes <- monocle3:::branch_nodes(ctl.cds, "UMAP")
branch_point_df <- ica_space_df %>%
  dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>%
  dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))
saveRDS(edge_df, file = paste0(rnaProject, "-ctl-learned_graph_DF.RDS"))
saveRDS(branch_point_df, file = paste0(rnaProject, "-ctl-learned_graph_branchpts.RDS"))

# Full object - T2D -------------------------------------------------------

t2d.down <- as.cell_data_set(t2d.down)
t2d.down <- cluster_cells(t2d.down)

t2d.down <- learn_graph(t2d.down)
saveRDS(t2d.down, file = paste0(rnaProject, "-t2d.down-t2d_learned-graph.RDS"))

# save graph components for plotting separately
ica_space_df <- t(t2d.down@principal_graph_aux$UMAP$dp_mst) %>%
  as.data.frame() %>%
  dplyr::select(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
  dplyr::mutate(sample_name = rownames(.),
                sample_state = rownames(.))
dp_mst <- t2d.down@principal_graph[["UMAP"]]
edge_df <- dp_mst %>%
  igraph::as_data_frame() %>%
  dplyr::select(source = "from", target = "to") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select(
                       source="sample_name",
                       source_prin_graph_dim_1="prin_graph_dim_1",
                       source_prin_graph_dim_2="prin_graph_dim_2"),
                   by = "source") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select(
                       target="sample_name",
                       target_prin_graph_dim_1="prin_graph_dim_1",
                       target_prin_graph_dim_2="prin_graph_dim_2"),
                   by = "target")

mst_branch_nodes <- monocle3:::branch_nodes(t2d.down, "UMAP")
branch_point_df <- ica_space_df %>%
  dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>%
  dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))
saveRDS(edge_df, file = paste0(rnaProject, "-t2d-learned_graph_DF.RDS"))
saveRDS(branch_point_df, file = paste0(rnaProject, "-t2d-learned_graph_branchpts.RDS"))

# Downsample split disease - Load object & subset  -----------------------

#load object
down.sample <- readRDS(file = paste0(rnaProject, "-downsampled.object.RDS"))
cds <- down.sample
cds[["UMAP"]] <- cds[["umap.rpca"]]
cds[["umap.rpca"]] <- NULL

#subset
ctl.down <- subset(cds, subset = Disease == "CTL")
t2d.down <- subset(cds, subset = Disease == "T2D")

# Downsample split disease - CTL  -----------------------

ctl.down <- as.cell_data_set(ctl.down)
ctl.down <- cluster_cells(ctl.down)

ctl.down <- learn_graph(ctl.down)
saveRDS(ctl.down, file = paste0(rnaProject, "-ctl.down-ctl_learned-graph.RDS"))

# save graph components for plotting separately
ica_space_df <- t(ctl.down@principal_graph_aux$UMAP$dp_mst) %>%
  as.data.frame() %>%
  dplyr::select(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
  dplyr::mutate(sample_name = rownames(.),
                sample_state = rownames(.))


dp_mst <- ctl.down@principal_graph[["UMAP"]]

edge_df <- dp_mst %>%
  igraph::as_data_frame() %>%
  dplyr::select(source = "from", target = "to") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select(
                       source="sample_name",
                       source_prin_graph_dim_1="prin_graph_dim_1",
                       source_prin_graph_dim_2="prin_graph_dim_2"),
                   by = "source") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select(
                       target="sample_name",
                       target_prin_graph_dim_1="prin_graph_dim_1",
                       target_prin_graph_dim_2="prin_graph_dim_2"),
                   by = "target")

mst_branch_nodes <- monocle3:::branch_nodes(ctl.down, "UMAP")
branch_point_df <- ica_space_df %>%
  dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>%
  dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))
saveRDS(edge_df, file = paste0(rnaProject, "-ctl-learned_graph_DF.RDS"))
saveRDS(branch_point_df, file = paste0(rnaProject, "-ctl-learned_graph_branchpts.RDS"))

# Downsample split Disease - T2D  -----------------------

t2d.down <- as.cell_data_set(t2d.down)
t2d.down <- cluster_cells(t2d.down)

t2d.down <- learn_graph(t2d.down)
saveRDS(t2d.down, file = paste0(rnaProject, "-t2d.down-t2d_learned-graph.RDS"))

# save graph components for plotting separately
ica_space_df <- t(t2d.down@principal_graph_aux$UMAP$dp_mst) %>%
  as.data.frame() %>%
  dplyr::select(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
  dplyr::mutate(sample_name = rownames(.),
                sample_state = rownames(.))
dp_mst <- t2d.down@principal_graph[["UMAP"]]
edge_df <- dp_mst %>%
  igraph::as_data_frame() %>%
  dplyr::select(source = "from", target = "to") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select(
                       source="sample_name",
                       source_prin_graph_dim_1="prin_graph_dim_1",
                       source_prin_graph_dim_2="prin_graph_dim_2"),
                   by = "source") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select(
                       target="sample_name",
                       target_prin_graph_dim_1="prin_graph_dim_1",
                       target_prin_graph_dim_2="prin_graph_dim_2"),
                   by = "target")

mst_branch_nodes <- monocle3:::branch_nodes(t2d.down, "UMAP")
branch_point_df <- ica_space_df %>%
  dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>%
  dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))
saveRDS(edge_df, file = paste0(rnaProject, "-t2d-learned_graph_DF.RDS"))
saveRDS(branch_point_df, file = paste0(rnaProject, "-t2d-learned_graph_branchpts.RDS"))
