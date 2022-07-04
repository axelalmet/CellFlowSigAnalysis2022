library(nichenetr)
library(tidyverse)
library(reticulate)
library(DiagrammeR)
library(igraph)
library(dplyr)
library(reticulate)

weighted_networks = readRDS("/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/NicheNet/nichenet_weighted_networks.rds")

ligand_tf_matrix = readRDS("/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/NicheNet/nichenet_ligand_tf_matrix.rds")

# Load the ligand-ligand edges
inferred_edges_liao <- read.csv("/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Liao2020/CausalityResults/liao20_cccflow_utigsp_adjacency_parcorr_bagged_ligand_target_edges.csv")

# We need to get the liao genes as well
ad <- import("anndata", convert = FALSE) # Import the module used to read in H5AD files
liao_ad_object <- ad$read_h5ad("/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Liao2020/liao20_sub.h5ad")
liao_genes <- rownames(py_to_r(liao_ad_object$var))

ligands_all = unique(inferred_edges_liao$Ligand1) # this can be a list of multiple ligands if required
targets_all = unique(inferred_edges_liao$Ligand2)

ligands_filtered = Reduce(intersect, list(rownames(ligand_tf_matrix), colnames(ligand_tf_matrix), ligands_all))
targets_filtered = intersect(union(weighted_networks_filtered$gr$from, weighted_networks_filtered$gr$to), targets_all)

# Filter the networks based on the genes present in liao etal.
ligand_tf_matrix_filtered <- ligand_tf_matrix[rownames(ligand_tf_matrix) %in% liao_genes, colnames(ligand_tf_matrix) %in% liao_genes]
weighted_networks_filtered <- weighted_networks
weighted_networks_filtered$lr_sig <- weighted_networks_filtered$lr_sig[(weighted_networks_filtered$lr_sig$from %in% liao_genes)|(weighted_networks_filtered$lr_sig$to %in% liao_genes), ]
weighted_networks_filtered$gr <- weighted_networks_filtered$gr[(weighted_networks_filtered$gr$from %in% liao_genes)|(weighted_networks_filtered$gr$to %in% liao_genes), ]

active_signaling_network_liao = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix_filtered, ligands_all = ligands_filtered, targets_all = targets_filtered, weighted_networks = weighted_networks_filtered)

# Construct network from active signaling network
liao_ligand_target_net <- graph_from_data_frame(d=rbind(active_signaling_network_liao$gr, active_signaling_network_liao$sig), directed=T) 

liao_vertices <- names(as.factor(V(liao_ligand_target_net)))
accepted_indices_liao <- c()
accepted_paths_liao <- c()
conditions_found_liao <- c()
for (i in 1:nrow(inferred_edges_liao))
{
  cluster_A <-  inferred_edges_liao[i, "Source"]
  ligand_A <- inferred_edges_liao[i, "Ligand1"]
  ligand_B <- inferred_edges_liao[i, "Ligand2"]
  receptor_B <- inferred_edges_liao[i, "Receptor"]
  cluster_B <-  inferred_edges_liao[i, "Target"]
  split_receptor_B <- strsplit(receptor_B, "_")[[1]]
  
  if ( (ligand_A %in% liao_vertices)&(ligand_B %in% liao_vertices) )
  {
    shortest_path_length <- 1e4
    shortest_path <- NULL
    
    for (j in 1:length(split_receptor_B))
    {
      sub <- split_receptor_B[j]
      if (sub %in% liao_vertices )
      {
        if (!is.infinite(distances(liao_ligand_target_net, sub, ligand_B, mode = "out")) )
        {
          current_path <- names(as.factor(all_shortest_paths(liao_ligand_target_net, sub, ligand_B, mode = "out")$res[[1]]))
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
      accepted_indices_liao <- c(accepted_indices_liao, i)
      
      # Get the shortest path from ligand to ligand via receptor and TFs
      accepted_paths_liao <- c(accepted_paths_liao, paste(c(ligand_A, shortest_path), collapse='_'))
      
    }
  }
}

implicated_edges_liao <- inferred_edges_liao[accepted_indices_liao, ]
implicated_edges_liao$Path <- accepted_paths_liao

write.csv(implicated_edges_liao, file = "/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Liao2020/CausalityResults/liao20_cccflow_utigsp_bagged_implicated_edges_by_nichenet.csv")
