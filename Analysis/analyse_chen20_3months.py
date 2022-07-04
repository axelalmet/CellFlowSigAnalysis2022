import scanpy as sc
import cellflowsig as cfs
import pandas as pd

data_directory = '/Users/axelalmet/Documents/CellCellCommunication/Data/Chen2020/'
cellchat_directory = '/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/CellChat/Output/Chen2020/'
giotto_directory = '/Users/axelalmet/Documents/CellCellCommunication/Data/Chen2020/'

# Load the scanpy object 
adata_chen = sc.read(data_directory + 'chen20_merged.h5ad')
adata_chen = adata_chen[adata_chen.obs['Group'].isin(['WT_03', 'AD_03'])]

# Import the cellchat dataframes
cellchat_chen_wt = pd.read_csv(cellchat_directory + 'chen20_communications_WT_03.csv')
cellchat_chen_ad = pd.read_csv(cellchat_directory + 'chen20_communications_AD_03.csv')

# Filter out weaker interactions to reduce noise.
# This step is optional, so if you don't want to do it, just comment out
# the four lines below
pval_threshold = 0.05
quantile = 0.25
cellchat_chen_wt = cellchat_chen_wt[(cellchat_chen_wt['pval'] < pval_threshold)&(cellchat_chen_wt['prob'] > cellchat_chen_wt['prob'].quantile(quantile))]
cellchat_chen_ad = cellchat_chen_ad[(cellchat_chen_ad['pval'] < pval_threshold)&(cellchat_chen_ad['prob'] > cellchat_chen_ad['prob'].quantile(quantile))]

# Define the cellchat output as a dictionary
cellchat_output_chen = {'WT_03':cellchat_chen_wt, 'AD_03':cellchat_chen_ad}

# Construct the list of feasible pairs based on spatial distances
# We manually encode the possible adjacencies (could probably do this from SquidPy, but i'm not sure how good that works)
enrichment_threshold = 0.0
feasible_pairs_chen = []

cellproximities_chen_ctrl = pd.read_csv(giotto_directory + 'chen20_cellproximities_coords_WT_03.csv')
cellproximities_chen_stim = pd.read_csv(giotto_directory + 'chen20_cellproximities_coords_AD_03.csv')
cellproximities_concat_chen = pd.concat([cellproximities_chen_ctrl, cellproximities_chen_stim])

# Any pair above the enrichment threshold is counted as a feasible pair
cellproximities_concat_chen = cellproximities_concat_chen[cellproximities_concat_chen['enrichm'] > enrichment_threshold]

for index, row in cellproximities_concat_chen.iterrows():

    cluster_A = row['cell_1']
    cluster_B = row['cell_2']

    # Add the pair, accounting for the fact that the distances are symmetric
    if ((cluster_A, cluster_B) not in feasible_pairs_chen)\
        |((cluster_A, cluster_B) not in feasible_pairs_chen):
        feasible_pairs_chen.append((cluster_A, cluster_B))

# Construct the base network now
cfs.pp.construct_base_networks(adata_chen,
                                cellchat_output_chen,
                                condition_label='Group',
                                method='cellchat',
                                feasible_pairs=feasible_pairs_chen,
                                node_sep='_',
                                base_network_label='base_networks')

# Construct the cell type ligand expression
cfs.pp.construct_celltype_ligand_expression(adata_chen,
                                            celltype_label='AT',
                                            node_sep='_',
                                            expressions_label='X_celltype_ligand',
                                            base_network_label='base_networks')

# Learn the causal signaling networks now
cfs.tl.learn_causal_network(adata_chen,
                            condition_label='Group',
                            control_label='WT_03',
                            causal_network_label='causal_networks',
                            celltype_ligands_label='X_celltype_ligand',
                            base_network_label='base_networks', 
                            n_cores=4,
                            n_bootstraps=100,
                            n_shuffles=1000,
                            alpha_ci=0.001,
                            alpha_inv=0.001)

# Validate causal network against base network
cfs.tl.validate_against_base_network(adata_chen,
                            condition_label='Group',
                            control_label='WT_03',
                            causal_network_label='causal_networks',
                            celltype_ligands_label='X_celltype_ligand',
                            base_network_label='base_networks')