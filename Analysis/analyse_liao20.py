import scanpy as sc
import cellflowsig as cfs
import pandas as pd

data_directory = '/Users/axelalmet/Documents/CellCellCommunication/Data/Liao2020/'
cellchat_directory = '/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/CellChat/Output/Liao2020/'
# Load teh 
adata_liao = sc.read(data_directory + 'liao20_merged.h5ad')


# Import the cellchat dataframes
cellchat_liao_healthy = pd.read_csv(cellchat_directory + 'liao20_communications_HC.csv')
cellchat_liao_moderate = pd.read_csv(cellchat_directory + 'liao20_communications_M.csv')
cellchat_liao_severe = pd.read_csv(cellchat_directory + 'liao20_communications_S.csv')

# Filter out weaker interactions to reduce noise.
# This step is optional, so if you don't want to do it, just comment out
# the four lines below
pval_threshold = 0.05
quantile = 0.25
cellchat_liao_healthy = cellchat_liao_healthy[(cellchat_liao_healthy['pval'] < pval_threshold)&(cellchat_liao_healthy['prob'] > cellchat_liao_healthy['prob'].quantile(quantile))]
cellchat_liao_moderate = cellchat_liao_moderate[(cellchat_liao_moderate['pval'] < pval_threshold)&(cellchat_liao_moderate['prob'] > cellchat_liao_moderate['prob'].quantile(quantile))]
cellchat_liao_severe = cellchat_liao_severe[(cellchat_liao_severe['pval'] < pval_threshold)&(cellchat_liao_severe['prob'] > cellchat_liao_severe['prob'].quantile(quantile))]

# Define the cellchat output as a dictionary
cellchat_output_liao = {'HC': cellchat_liao_healthy, 'M':cellchat_liao_moderate, 'S':cellchat_liao_severe}

# Construct the base network now
cfs.pp.construct_base_networks(adata_liao,
                                cellchat_output_liao,
                                condition_label='Condition',
                                method='cellchat',
                                node_sep='_',
                                base_network_label='base_networks')

# Construct the cell type ligand expression
cfs.pp.construct_celltype_ligand_expressions(adata_liao,
                                            celltype_label='celltype',
                                            node_sep='_',
                                            expressions_label='X_celltype_ligand',
                                            base_network_label='base_networks')

# Learn the causal signaling networks now
cfs.tl.learn_causal_network(adata_liao,
                            condition_label='group',
                            control_label='HC',
                            causal_network_label='causal_networks',
                            celltype_ligands_label='X_celltype_ligand',
                            base_network_label='base_networks', 
                            n_cores=8,
                            n_bootstraps=100,
                            n_shuffles=1000,
                            alpha_ci=0.001,
                            alpha_inv=0.001)

# Validate causal network against base network
cfs.tl.validate_against_base_network(adata_liao,
                            condition_label='Condition',
                            causal_network_label='causal_networks',
                            celltype_ligands_label='X_celltype_ligand',
                            base_network_label='base_networks')