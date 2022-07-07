import scanpy as sc
import cellflowsig as cfs
import pandas as pd

data_directory = '/Users/axelalmet/Documents/CellFlowSigAnalysis2022/SyntheticBenchmark/Fork/'
cellchat_directory = data_directory + 'Output/CellCellCommunication/'

# Load the scRNA-seq data
adata_fork = sc.read(data_directory + 'adata_fork.h5ad')

# Import the cellchat dataframes
cellchat_fork_obs = pd.read_csv(cellchat_directory + 'fork_communications_observation.csv')
cellchat_fork_int = pd.read_csv(cellchat_directory + 'fork_communications_intervention.csv')

# Filter out weaker interactions to reduce noise.
# This step is optional, so if you don't want to do it, just comment out
# the four lines below
pval_threshold = 0.05
quantile = 0.25
cellchat_fork_obs = cellchat_fork_obs[(cellchat_fork_obs['pval'] < pval_threshold)&(cellchat_fork_obs['prob'] > cellchat_fork_obs['prob'].quantile(quantile))]
cellchat_fork_int = cellchat_fork_int[(cellchat_fork_int['pval'] < pval_threshold)&(cellchat_fork_int['prob'] > cellchat_fork_int['prob'].quantile(quantile))]

# Define the cellchat output as a dictionary
cellchat_output_fork = {'Observation':cellchat_fork_obs, 'Intervention':cellchat_fork_int}

# Construct the base network now
cfs.pp.construct_base_networks(adata_fork,
                                cellchat_output_fork,
                                condition_label='Condition',
                                method='cellchat',
                                node_sep='_',
                                base_network_label='base_networks')

# Construct the cell type ligand expression
cfs.pp.construct_celltype_ligand_expressions(adata_fork,
                                            celltype_label='type',
                                            node_sep='_',
                                            expressions_label='X_celltype_ligand',
                                            base_network_label='base_networks')

# Learn the causal signaling networks now
cfs.tl.learn_causal_network(adata_fork,
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

# print(adata_fork.uns['causal_networks'])

# # Validate causal network against base network
cfs.tl.validate_against_base_network(adata_fork,
                                    condition_label = 'Condition',
                                    causal_network_label = 'causal_networks',
                                    base_network_label  = 'base_networks')

print(adata_fork.uns['causal_networks'])
