import scanpy as sc
import cellflowsig as cfs
import pandas as pd

data_directory = '/Users/axelalmet/Documents/CellCellCommunication/Data/Mcfaline-Figueroa2019/'
cellchat_directory = '/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/CellChat/Output/Mcfaline-Figueroa2019/'
# Load teh 
adata_mcfaline = sc.read(data_directory + 'mcfaline19_merged.h5ad')

# Import the cellchat dataframes
cellchat_mcfaline_inner = pd.read_csv(cellchat_directory + 'mcfaline19_communications_inner.csv')
cellchat_mcfaline_outer = pd.read_csv(cellchat_directory + 'mcfaline19_communications_outer.csv')

# Filter out weaker interactions to reduce noise.
# This step is optional, so if you don't want to do it, just comment out
# the four lines below
pval_threshold = 0.05
quantile = 0.25
cellchat_mcfaline_inner = cellchat_mcfaline_inner[(cellchat_mcfaline_inner['pval'] < pval_threshold)&(cellchat_mcfaline_inner['prob'] > cellchat_mcfaline_inner['prob'].quantile(quantile))]
cellchat_mcfaline_outer = cellchat_mcfaline_outer[(cellchat_mcfaline_outer['pval'] < pval_threshold)&(cellchat_mcfaline_outer['prob'] > cellchat_mcfaline_outer['prob'].quantile(quantile))]

# Define the cellchat output as a list
cellchat_output_mcfaline = {'inner':cellchat_mcfaline_inner, 'outer':cellchat_mcfaline_outer}

# Construct the base network now
cfs.pp.construct_base_networks(adata_mcfaline,
                                cellchat_output_mcfaline,
                                condition_label='spatial_id',
                                method='cellchat',
                                celltype_sep_old=' ',
                                celltype_sep_new='-', 
                                node_sep='_',
                                base_network_label='base_networks')

# Construct the cell type ligand expression
cfs.pp.construct_celltype_ligand_expressions(adata_mcfaline,
                                            celltype_label='leiden',
                                            celltype_sep_old=' ',
                                            celltype_sep_new='-',
                                            node_sep='_',
                                            expressions_label='X_celltype_ligand',
                                            base_network_label='base_networks')

# Learn the causal signaling networks now
cfs.tl.learn_causal_network(adata_mcfaline,
                            condition_label='spatial_id',
                            control_label='inner',
                            causal_network_label='causal_networks',
                            celltype_ligands_label='X_celltype_ligand',
                            base_network_label='base_networks', 
                            n_cores=4,
                            n_bootstraps=100,
                            n_shuffles=1000,
                            alpha_ci=0.001,
                            alpha_inv=0.001)

# Validate causal network against base network
cfs.tl.validate_against_base_network(adata_mcfaline,
                            condition_label='spatial_id',
                            causal_network_label='causal_networks',
                            celltype_ligands_label='X_celltype_ligand',
                            base_network_label='base_networks',
                            celltype_sep_old=' ',
                            celltype_sep_new='-',)