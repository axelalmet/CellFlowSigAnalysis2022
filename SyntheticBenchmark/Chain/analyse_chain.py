import scanpy as sc
import cellflowsig as cfs
import pandas as pd

data_directory = '/Users/axelalmet/Documents/CellFlowSigAnalysis2022/SyntheticBenchmark/Chain/'
cellchat_directory = data_directory + 'Output/CellCellCommunication/'

# Load the scRNA-seq data
adata_chain = sc.read(data_directory + 'adata_chain.h5ad')

# Import the cellchat dataframes
cellchat_chain_obs = pd.read_csv(cellchat_directory + 'chain_communications_observation.csv')
cellchat_chain_int = pd.read_csv(cellchat_directory + 'chain_communications_intervention.csv')

# Filter out weaker interactions to reduce noise.
# This step is optional, so if you don't want to do it, just comment out
# the four lines below
pval_threshold = 0.05
quantile = 0.25
cellchat_chain_obs = cellchat_chain_obs[(cellchat_chain_obs['pval'] < pval_threshold)&(cellchat_chain_obs['prob'] > cellchat_chain_obs['prob'].quantile(quantile))]
cellchat_chain_int = cellchat_chain_int[(cellchat_chain_int['pval'] < pval_threshold)&(cellchat_chain_int['prob'] > cellchat_chain_int['prob'].quantile(quantile))]

# Define the cellchat output as a dictionary
cellchat_output_chain = {'Observation':cellchat_chain_obs, 'Intervention':cellchat_chain_int}

# Construct the base network now
cfs.pp.construct_base_networks(adata_chain,
                                cellchat_output_chain,
                                condition_label='Condition',
                                method='cellchat',
                                node_sep='_',
                                base_network_label='base_networks')

# Construct the cell type ligand expression
cfs.pp.construct_celltype_ligand_expressions(adata_chain,
                                            celltype_label='type',
                                            node_sep='_',
                                            expressions_label='X_celltype_ligand',
                                            base_network_label='base_networks')

# Learn the causal signaling networks now
cfs.tl.learn_causal_network(adata_chain,
                            condition_label='Condition',
                            control_label='Observation',
                            causal_network_label='causal_networks',
                            celltype_ligands_label='X_celltype_ligand',
                            base_network_label='base_networks', 
                            n_jobs=4,
                            n_bootstraps=100,
                            n_shuffles=1000,
                            alpha_ci=0.001,
                            alpha_inv=0.001)

# print(adata_chain.uns['causal_networks'])

# # Validate causal network against base network
cfs.tl.validate_against_base_network(adata_chain,
                                    condition_label = 'Condition',
                                    causal_network_label = 'causal_networks',
                                    base_network_label  = 'base_networks')

print(adata_chain.uns['causal_networks'])
