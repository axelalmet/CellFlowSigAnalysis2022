library(nichenetr)
library(tidyverse)
library(reticulate)
library(DiagrammeR)
library(igraph)
library(dplyr)
library(reticulate)

weighted_networks = readRDS("/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/NicheNet/nichenet_weighted_networks.rds")

weighted_networks$lr_sig$from <- convert_human_to_mouse_symbols(weighted_networks$lr_sig$from)
weighted_networks$lr_sig$to <- convert_human_to_mouse_symbols(weighted_networks$lr_sig$to)
weighted_networks$lr_sig <- weighted_networks$lr_sig[(!is.na(weighted_networks$lr_sig$from)&!is.na(weighted_networks$lr_sig$to)), ]

weighted_networks$gr$from <- convert_human_to_mouse_symbols(weighted_networks$gr$from)
weighted_networks$gr$to <- convert_human_to_mouse_symbols(weighted_networks$gr$to)
weighted_networks$gr <- weighted_networks$gr[(!is.na(weighted_networks$gr$from)&!is.na(weighted_networks$gr$to)), ]

ligand_tf_matrix = readRDS("/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/NicheNet/nichenet_ligand_tf_matrix.rds")

colnames(ligand_tf_matrix) = ligand_tf_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_tf_matrix) = ligand_tf_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
ligand_tf_matrix = ligand_tf_matrix %>% .[!is.na(rownames(ligand_tf_matrix)), !is.na(colnames(ligand_tf_matrix))]

# Load the ligand-ligand edges
inferred_edges_chen <- read.csv("/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Chen2020/CausalityResults/chen20_3mo_cccflow_utigsp_adjacency_parcorr_bagged_ligand_target_edges.csv")

# We need to get the chen genes as well
ad <- import("anndata", convert = FALSE) # Import the module used to read in H5AD files
chen_ad_object <- ad$read_h5ad("/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Chen2020/chen20_merged.h5ad")
chen_genes <- rownames(py_to_r(chen_ad_object$var))

ligands_all = unique(inferred_edges_chen$Ligand1) # this can be a list of multiple ligands if required
targets_all = unique(inferred_edges_chen$Ligand2)

ligands_filtered = Reduce(intersect, list(rownames(ligand_tf_matrix), colnames(ligand_tf_matrix), ligands_all))
targets_filtered = intersect(union(weighted_networks$gr$from, weighted_networks$gr$to), targets_all)

# Filter the networks based on the genes present in chen etal.
ligand_tf_matrix_filtered <- ligand_tf_matrix[rownames(ligand_tf_matrix) %in% chen_genes, colnames(ligand_tf_matrix) %in% chen_genes]
weighted_networks_filtered <- weighted_networks
weighted_networks_filtered$lr_sig <- weighted_networks_filtered$lr_sig[(weighted_networks_filtered$lr_sig$from %in% chen_genes)|(weighted_networks_filtered$lr_sig$to %in% chen_genes), ]
weighted_networks_filtered$gr <- weighted_networks_filtered$gr[(weighted_networks_filtered$gr$from %in% chen_genes)|(weighted_networks_filtered$gr$to %in% chen_genes), ]

active_signaling_network_chen = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix_filtered, ligands_all = ligands_filtered, targets_all = targets_filtered, weighted_networks = weighted_networks_filtered)
# Construct network from active signaling network
chen_ligand_target_net <- graph_from_data_frame(d=rbind(active_signaling_network_chen$gr, active_signaling_network_chen$sig), directed=T) 

chen_vertices <- names(as.factor(V(chen_ligand_target_net)))
accepted_indices_chen <- c()
accepted_paths_chen <- c()
for (i in 1:nrow(inferred_edges_chen))
{
  cluster_A <- inferred_edges_chen[i, "Source"]
  ligand_A <- inferred_edges_chen[i, "Ligand1"]
  ligand_B <- inferred_edges_chen[i, "Ligand2"]
  receptor_B <- inferred_edges_chen[i, "Receptor"]
  cluster_B <- inferred_edges_chen[i, "Target"]
  split_receptor_B <- strsplit(receptor_B, "_")[[1]]
  
  shortest_path_length <- 1e4
  shortest_path <- NULL
  for (j in 1:length(split_receptor_B))
  {
    sub <- split_receptor_B[j]
    if (sub %in% chen_vertices )
    {
      if (!is.infinite(distances(chen_ligand_target_net, sub, ligand_B, mode = "out")) )
      {
        current_path <- names(as.factor(all_shortest_paths(chen_ligand_target_net, sub, ligand_B, mode = "out")$res[[1]]))
        if (length(current_path) < shortest_path_length)
        {
          shortest_path_length <- length(current_path)
          shortest_path <- paste(current_path, collapse='_')
        }
      }
    }
  }
  if (!is.null(shortest_path)) # If the shortest path length changed, that means we can accept the index
  {
    accepted_indices_chen <- c(accepted_indices_chen, i)
    
    # Get the shortest path from ligand to ligand via receptor and TFs
    accepted_paths_chen <- c(accepted_paths_chen, paste(c(ligand_A, shortest_path), collapse='_'))
    
  }
  
}

implicated_edges_chen <- inferred_edges_chen[accepted_indices_chen, ]
implicated_edges_chen$Path <- accepted_paths_chen

write.csv(implicated_edges_chen, file = "/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Chen2020/CausalityResults/chen20_3mo_cccflow_utigsp_bagged_implicated_edges_by_nichenet.csv")
