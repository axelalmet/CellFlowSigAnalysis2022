library(nichenetr)
library(tidyverse)
library(reticulate)
library(DiagrammeR)
library(igraph)
library(dplyr)

weighted_networks = readRDS("/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/NicheNet/nichenet_weighted_networks.rds")

ligand_tf_matrix = readRDS("/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/NicheNet/nichenet_ligand_tf_matrix.rds")

# Load the ligan
inferred_edges_kang <- read.csv("/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Kang2018/CausalityResults/kang18_cccflow_utigsp_adjacency_parcorr_bagged_ligand_target_edges.csv")

# We need to get the kang genes as well
ad <- import("anndata", convert = FALSE) # Import the module used to read in H5AD files
kang_ad_object <- ad$read_h5ad("/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Kang2018/kang18_tutorial.h5ad")
kang_genes <- rownames(py_to_r(kang_ad_object$var))

ligands_all = unique(inferred_edges_kang$Ligand1) # this can be a list of multiple ligands if required
targets_all = unique(inferred_edges_kang$Ligand2)

ligands_filtered = Reduce(intersect, list(rownames(ligand_tf_matrix), colnames(ligand_tf_matrix), ligands_all))

# Filter the networks based on the genes present in kang etal.
ligand_tf_matrix_filtered <- ligand_tf_matrix[rownames(ligand_tf_matrix) %in% kang_genes, colnames(ligand_tf_matrix) %in% kang_genes]
weighted_networks_filtered <- weighted_networks
weighted_networks_filtered$lr_sig <- weighted_networks_filtered$lr_sig[(weighted_networks_filtered$lr_sig$from %in% kang_genes)|(weighted_networks_filtered$lr_sig$to %in% kang_genes), ]
weighted_networks_filtered$gr <- weighted_networks_filtered$gr[(weighted_networks_filtered$gr$from %in% kang_genes)|(weighted_networks_filtered$gr$to %in% kang_genes), ]

active_signaling_network_kang = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix_filtered, ligands_all = ligands_filtered, targets_all = targets_all, weighted_networks = weighted_networks_filtered)

# Construct network from active signaling network
kang_ligand_target_net <- graph_from_data_frame(d=rbind(active_signaling_network_kang$gr, active_signaling_network_kang$sig), directed=T) 

kang_vertices <- names(as.factor(V(kang_ligand_target_net)))
accepted_indices_kang <- c()
accepted_paths_kang <- c()
for (i in 1:nrow(inferred_edges_kang))
{
  cluster_A <- inferred_edges_kang[i, "Source"]
  ligand_A <- inferred_edges_kang[i, "Ligand1"]
  ligand_B <- inferred_edges_kang[i, "Ligand2"]
  receptor_B <- inferred_edges_kang[i, "Receptor"]
  cluster_B <-  inferred_edges_kang[i, "Target"]
  split_receptor_B <- strsplit(receptor_B, "_")[[1]]
  
  for (j in 1:length(split_receptor_B))
  {
    sub <- split_receptor_B[j]
    if (sub %in% kang_vertices )
    {
      if (!is.infinite(distances(kang_ligand_target_net, sub, ligand_B, mode = "out")) )
      {
        if (!i %in% accepted_indices_kang)
        {
          accepted_indices_kang <- c(accepted_indices_kang, i)
          
          # Get the shortest path from ligand to ligand via receptor and TFs
          shortest_path <- paste(names(as.factor(all_shortest_paths(kang_ligand_target_net, sub, ligand_B, mode = "out")$res[[1]])), collapse='_')
          accepted_paths_kang <- c(accepted_paths_kang, paste(c(ligand_A, shortest_path), collapse='_'))
          
        }
      }
      
    }
  }
  
}


implicated_edges_kang <- inferred_edges_kang[accepted_indices_kang, ]
implicated_edges_kang$Path <- accepted_paths_kang

write.csv(implicated_edges_kang, file = "/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Kang2018/CausalityResults/kang18_cccflow_utigsp_bagged_implicated_edges_by_nichenet.csv")
