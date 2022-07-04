library(nichenetr)
library(tidyverse)
library(reticulate)
library(igraph)
library(dplyr)
library(reticulate)

# Load all the necessary matrices and convert the human genes to mouse symbols
weighted_networks = readRDS("../Networks/nichenet_weighted_networks.rds")

ligand_tf_matrix = readRDS("/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/NicheNet/nichenet_ligand_tf_matrix.rds")

# Load the ligand-target edges
inferred_edges_mcfaline <- read.csv("/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/McFaline-Figueroa2019/CausalityResults/mcfaline19_cccflow_utigsp_adjacency_parcorr_bagged_ligand_target_edges.csv")

# We need to get the mcfaline genes as well
ad <- import("anndata", convert = FALSE) # Import the module used to read in H5AD files
mcfaline_ad_object <- ad$read_h5ad("/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/McFaline-Figueroa2019/mcfaline2019_merged.h5ad")
mcfaline_genes <- rownames(py_to_r(mcfaline_ad_object$var))

ligands_all = unique(inferred_edges_mcfaline$Ligand1) # this can be a list of multiple ligands if required
targets_all = unique(inferred_edges_mcfaline$Ligand2)

ligands_filtered = Reduce(intersect, list(rownames(ligand_tf_matrix), colnames(ligand_tf_matrix), ligands_all))

# Filter the networks based on the genes present in mcfaline etal.
ligand_tf_matrix_filtered <- ligand_tf_matrix[rownames(ligand_tf_matrix) %in% mcfaline_genes, colnames(ligand_tf_matrix) %in% mcfaline_genes]
weighted_networks_filtered <- weighted_networks
weighted_networks_filtered$lr_sig <- weighted_networks_filtered$lr_sig[(weighted_networks_filtered$lr_sig$from %in% mcfaline_genes)|(weighted_networks_filtered$lr_sig$to %in% mcfaline_genes), ]
weighted_networks_filtered$gr <- weighted_networks_filtered$gr[(weighted_networks_filtered$gr$from %in% mcfaline_genes)|(weighted_networks_filtered$gr$to %in% mcfaline_genes), ]

active_signaling_network_mcfaline = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix_filtered, ligands_all = ligands_filtered, targets_all = targets_all, weighted_networks = weighted_networks_filtered)

# Construct network from active signaling network
mcfaline_ligand_target_net <- graph_from_data_frame(d=rbind(active_signaling_network_mcfaline$gr, active_signaling_network_mcfaline$sig), directed=T) 

mcfaline_vertices <- names(as.factor(V(mcfaline_ligand_target_net)))
accepted_indices_mcfaline <- c()
accepted_paths_mcfaline <- c()
for (i in 1:nrow(inferred_edges_mcfaline))
{
  cluster_A <-  inferred_edges_mcfaline[i, "Source"]
  ligand_A <- inferred_edges_mcfaline[i, "Ligand1"]
  ligand_B <- inferred_edges_mcfaline[i, "Ligand2"]
  receptor_B <- inferred_edges_mcfaline[i, "Receptor"]
  cluster_B <-  inferred_edges_mcfaline[i, "Target"]
  split_receptor_B <- strsplit(receptor_B, "_")[[1]]
  
  shortest_path_length <- 1e4
  shortest_path <- NULL
  for (j in 1:length(split_receptor_B))
  {
    sub <- split_receptor_B[j]
    if (sub %in% mcfaline_vertices )
    {
      if (!is.infinite(distances(mcfaline_ligand_target_net, sub, ligand_B, mode = "out")) )
      {
        current_path <- names(as.factor(all_shortest_paths(mcfaline_ligand_target_net, sub, ligand_B, mode = "out")$res[[1]]))
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
    accepted_indices_mcfaline <- c(accepted_indices_mcfaline, i)
    
    # Get the shortest path from ligand to ligand via receptor and TFs
    accepted_paths_mcfaline <- c(accepted_paths_mcfaline, paste(c(ligand_A, shortest_path), collapse='_'))

  }
  
}

implicated_edges_mcfaline <- inferred_edges_mcfaline[accepted_indices_mcfaline, ]
implicated_edges_mcfaline$Path <- accepted_paths_mcfaline
write.csv(implicated_edges_mcfaline, file = "/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/McFaline-Figueroa2019/CausalityResults/mcfaline19_cccflow_utigsp_bagged_implicated_edges_by_nichenet.csv")
