import scanpy as sc
import cellflowsig as cfs
import pandas as pd

data_directory = '/Users/axelalmet/Documents/CellCellCommunication/Data/Burkhardt2021/'
cellchat_directory = '/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/CellChat/Output/Burkhardt2021/'
# Load teh 
adata_burkhardt = sc.read(data_directory + 'burkhardt21_merged.h5ad')


# Import the cellchat dataframes
cellchat_burkhardt_ctrl = pd.read_csv(cellchat_directory + 'burkhardt21_communications_Ctrl.csv')
cellchat_burkhardt_ifng = pd.read_csv(cellchat_directory + 'burkhardt21_communications_IFNg.csv')

# Filter out weaker interactions to reduce noise.
# This step is optional, so if you don't want to do it, just comment out
# the four lines below
pval_threshold = 0.05
quantile = 0.25
cellchat_burkhardt_ctrl = cellchat_burkhardt_ctrl[(cellchat_burkhardt_ctrl['pval'] < pval_threshold)&(cellchat_burkhardt_ctrl['prob'] > cellchat_burkhardt_ctrl['prob'].quantile(quantile))]
cellchat_burkhardt_ifng = cellchat_burkhardt_ifng[(cellchat_burkhardt_ifng['pval'] < pval_threshold)&(cellchat_burkhardt_ifng['prob'] > cellchat_burkhardt_ifng['prob'].quantile(quantile))]

# Define the cellchat output as a dictionary
cellchat_output_burkhardt = {'Ctrl':cellchat_burkhardt_ctrl, 'IFNg':cellchat_burkhardt_ifng}

# Construct the base network now
cfs.pp.construct_base_networks(adata_burkhardt,
                                cellchat_output_burkhardt,
                                condition_label='Condition',
                                method='cellchat',
                                node_sep='_',
                                base_network_label='base_networks')

# Construct the cell type ligand expression
cfs.pp.construct_celltype_ligand_expressions(adata_burkhardt,
                                            celltype_label='Type',
                                            node_sep='_',
                                            expressions_label='X_celltype_ligand',
                                            base_network_label='base_networks')

# Learn the causal signaling networks now
cfs.tl.learn_causal_network(adata_burkhardt,
                            condition_label='Condition',
                            control_label='Ctrl',
                            causal_network_label='causal_networks',
                            celltype_ligands_label='X_celltype_ligand',
                            base_network_label='base_networks', 
                            n_cores=4,
                            n_bootstraps=100,
                            n_shuffles=1000,
                            alpha_ci=0.001,
                            alpha_inv=0.001)

# Validate causal network against base network
cfs.tl.validate_against_base_network(adata_burkhardt,
                            condition_label='Condition',
                            control_label='Ctrl',
                            causal_network_label='causal_networks',
                            celltype_ligands_label='X_celltype_ligand',
                            base_network_label='base_networks')