library(nichenetr)
library(tidyverse)
library(reticulate)
library(DiagrammeR)
library(igraph)
library(dplyr)

weighted_networks = readRDS("/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/NicheNet/nichenet_weighted_networks.rds")
# weighted_networks$lr_sig$from <- convert_human_to_mouse_symbols(weighted_networks$lr_sig$from)
# weighted_networks$lr_sig$to <- convert_human_to_mouse_symbols(weighted_networks$lr_sig$to)
# weighted_networks$lr_sig <- weighted_networks$lr_sig[(!is.na(weighted_networks$lr_sig$from)&!is.na(weighted_networks$lr_sig$to)), ]

# weighted_networks$gr$from <- convert_human_to_mouse_symbols(weighted_networks$gr$from)
# weighted_networks$gr$to <- convert_human_to_mouse_symbols(weighted_networks$gr$to)
# weighted_networks$gr <- weighted_networks$gr[(!is.na(weighted_networks$gr$from)&!is.na(weighted_networks$gr$to)), ]

ligand_tf_matrix = readRDS("/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/NicheNet/nichenet_ligand_tf_matrix.rds")
# colnames(ligand_tf_matrix) = ligand_tf_matrix %>% colnames() %>% convert_human_to_mouse_symbols() 
# rownames(ligand_tf_matrix) = ligand_tf_matrix %>% rownames() %>% convert_human_to_mouse_symbols() 
# ligand_tf_matrix = ligand_tf_matrix %>% .[!is.na(rownames(ligand_tf_matrix)), !is.na(colnames(ligand_tf_matrix))]

# Load the ligand-target relations
inferred_edges_burkhardt <- read.csv("/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Burkhardt2021/CausalityResults/burkhardt21_cccflow_utigsp_adjacency_parcorr_bagged_ligand_target_edges.csv")

# We need to get the burkhardt genes as well
ad <- import("anndata", convert = FALSE) # Import the module used to read in H5AD files
burkhardt_ad_object <- ad$read_h5ad("/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Burkhardt2021/burkhardt_merged.h5ad")
burkhardt_genes <- rownames(py_to_r(burkhardt_ad_object$var))

ligands_all = unique(inferred_edges_burkhardt$Ligand1) # this can be a list of multiple ligands if required
targets_all = unique(inferred_edges_burkhardt$Ligand2)

# Filter the networks based on the genes present in burkhardt etal.
ligand_tf_matrix_filtered <- ligand_tf_matrix[rownames(ligand_tf_matrix) %in% burkhardt_genes, colnames(ligand_tf_matrix) %in% burkhardt_genes]
weighted_networks_filtered <- weighted_networks
weighted_networks_filtered$lr_sig <- weighted_networks_filtered$lr_sig[(weighted_networks_filtered$lr_sig$from %in% burkhardt_genes)|(weighted_networks_filtered$lr_sig$to %in% burkhardt_genes), ]
weighted_networks_filtered$gr <- weighted_networks_filtered$gr[(weighted_networks_filtered$gr$from %in% burkhardt_genes)|(weighted_networks_filtered$gr$to %in% burkhardt_genes), ]

active_signaling_network_burkhardt = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix_filtered, ligands_all = ligands_all, targets_all = targets_all, weighted_networks = weighted_networks_filtered)

# Construct network from active signaling network
burkhardt_ligand_target_net <- graph_from_data_frame(d=rbind(active_signaling_network_burkhardt$sig, active_signaling_network_burkhardt$gr), directed=T) 

burkhardt_vertices <- names(as.factor(V(burkhardt_ligand_target_net)))
accepted_indices_burkhardt <- c()
accepted_paths_burkhardt <- c()
conditions_found_burkhardt <- c()
for (i in 1:nrow(inferred_edges_burkhardt))
{
  cluster_A <- inferred_edges_burkhardt[i, "Source"]
  ligand_A <- inferred_edges_burkhardt[i, "Ligand1"]
  ligand_B <- inferred_edges_burkhardt[i, "Ligand2"]
  receptor_B <- inferred_edges_burkhardt[i, "Receptor"]
  cluster_B <- inferred_edges_burkhardt[i, "Target"]
  split_receptor_B <- strsplit(receptor_B, "_")[[1]]

  shortest_path_length <- 1e4
  shortest_path <- NULL
  for (j in 1:length(split_receptor_B))
  {
    sub <- split_receptor_B[j]
    if (sub %in% burkhardt_vertices )
    {
      if (!is.infinite(distances(burkhardt_ligand_target_net, sub, ligand_B, mode = "out")) )
      {
        current_path <- names(as.factor(all_shortest_paths(burkhardt_ligand_target_net, sub, ligand_B, mode = "out")$res[[1]]))
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
    accepted_indices_burkhardt <- c(accepted_indices_burkhardt, i)
    
    # Get the shortest path from ligand to ligand via receptor and TFs
    accepted_paths_burkhardt <- c(accepted_paths_burkhardt, paste(c(ligand_A, shortest_path), collapse='_'))
    
  }
}

implicated_edges_burkhardt <- inferred_edges_burkhardt[accepted_indices_burkhardt, ]
implicated_edges_burkhardt$Path <- accepted_paths_burkhardt
write.csv(implicated_edges_burkhardt, file = "/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Burkhardt2021/CausalityResults/burkhardt21_cccflow_utigsp_bagged_implicated_edges_by_nichenet.csv")
