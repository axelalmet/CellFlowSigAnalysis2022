import scanpy as sc
import cellflowsig as cfs
import pandas as pd

data_directory = '/Users/axelalmet/Documents/CellCellCommunicationData/'
cellchat_directory = '/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/CellChat/Output/Kang2018/'
# Load teh 
adata_kang = sc.read(data_directory + 'kang18_merged.h5ad')


# Import the cellchat dataframes
cellchat_kang_ctrl = pd.read_csv(cellchat_directory + 'kang18_communications_control.csv')
cellchat_kang_stim = pd.read_csv(cellchat_directory + 'kang18_communications_stimulated.csv')

# Filter out weaker interactions to reduce noise.
# This step is optional, so if you don't want to do it, just comment out
# the four lines below
pval_threshold = 0.05
quantile = 0.25
cellchat_kang_ctrl = cellchat_kang_ctrl[(cellchat_kang_ctrl['pval'] < pval_threshold)&(cellchat_kang_ctrl['prob'] > cellchat_kang_ctrl['prob'].quantile(quantile))]
cellchat_kang_stim = cellchat_kang_stim[(cellchat_kang_stim['pval'] < pval_threshold)&(cellchat_kang_stim['prob'] > cellchat_kang_stim['prob'].quantile(quantile))]

# Define the cellchat output as a dictionary
cellchat_output_kang = {'control':cellchat_kang_ctrl, 'stimulated':cellchat_kang_stim}

# Construct the base network now
cfs.pp.construct_base_networks(adata_kang,
                                cellchat_output_kang,
                                condition_label='condition',
                                method='cellchat',
                                celltype_sep_old=' ',
                                celltype_sep_new='-', 
                                node_sep='_',
                                base_network_label='base_networks')

# Construct the cell type ligand expression
cfs.pp.construct_celltype_ligand_expressions(adata_kang,
                                            celltype_label='cell_type',
                                            celltype_sep_old=' ',
                                            celltype_sep_new='-', 
                                            node_sep='_',
                                            expressions_label='X_celltype_ligand',
                                            base_network_label='base_networks')

# Learn the causal signaling networks now
cfs.tl.learn_causal_network(adata_kang,
                            condition_label='condition',
                            control_label='control',
                            causal_network_label='causal_networks',
                            celltype_ligands_label='X_celltype_ligand',
                            base_network_label='base_networks', 
                            n_cores=4,
                            n_bootstraps=100,
                            n_shuffles=1000,
                            alpha_ci=0.001,
                            alpha_inv=0.001)

# Validate causal network against base network
cfs.tl.validate_against_base_network(adata_kang,
                            condition_label='condition',
                            control_label='control',
                            causal_network_label='causal_networks',
                            celltype_ligands_label='X_celltype_ligand',
                            base_network_label='base_networks')