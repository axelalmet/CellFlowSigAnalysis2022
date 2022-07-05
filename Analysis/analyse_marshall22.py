import scanpy as sc
import cellflowsig as cfs
import pandas as pd

data_directory = '/Users/axelalmet/Documents/CellCellCommunication/Data/Marshall2022/'
cellchat_directory = '/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/CellChat/Output/Marshall2022/'
giotto_directory = '/Users/axelalmet/Documents/CellCellCommunication/Data/Marshall2022/'

# Load the scanpy object 
adata_marshall = sc.read(data_directory + 'marshall22_umod_merged.h5ad')

# We construct the cellchat output and feasible pairs list simultaneously.
# To deal with the noise in the Slide-SeqV2 data, we construct the cellchat
# and cell proximity networks on a per-sub-sample basis.
enrichment_threshold = 0.0
feasible_pairs_marshall = []

WT_samples = sorted(adata_marshall[adata_marshall.obs['sample'] == 'WT'].obs['sub_sample'].unique().tolist())
KI_samples = sorted(adata_marshall[adata_marshall.obs['sample'] == 'KI'].obs['sub_sample'].unique().tolist())
cellproximities_marshall_wt_list = []
cellproximities_marshall_ki_list = []
cellchat_marshall_wt_list = []
cellchat_marshall_ki_list = []

# Import the cell proximities and cellchat objects for each list
for i in range(len(WT_samples)):
    wt_sub_sample = WT_samples[i]
    ki_sub_sample = KI_samples[i]

    cellproximities_marshall_wt = pd.read_csv(giotto_directory + 'marshall22_umod_cellproximities_coords_' + wt_sub_sample + '.csv')
    cellproximities_marshall_ki = pd.read_csv(giotto_directory + 'marshall22_umod_cellproximities_coords_' + ki_sub_sample + '.csv')
    cellchat_marshall_wt = pd.read_csv(cellchat_directory + 'marshall22_umod_communications_' + wt_sub_sample + '.csv')
    cellchat_marshall_ki = pd.read_csv(cellchat_directory + 'marshall22_umod_communications_' + ki_sub_sample + '.csv')

    cellproximities_marshall_wt_list.append(cellproximities_marshall_wt)
    cellproximities_marshall_ki_list.append(cellproximities_marshall_ki)
    cellchat_marshall_wt_list .append(cellchat_marshall_wt)
    cellchat_marshall_ki_list.append(cellchat_marshall_ki)

# Concatenate the cell proximities lists
cellproximities_marshall_wt = pd.concat(cellproximities_marshall_wt_list)
cellproximities_marshall_ki = pd.concat(cellproximities_marshall_ki_list)
cellproximities_concat_marshall = pd.concat([cellproximities_marshall_wt, cellproximities_marshall_ki])

# Filter out edges below the enrichment threshold
cellproximities_concat_marshall = cellproximities_concat_marshall[cellproximities_concat_marshall['enrichm'] > enrichment_threshold]

for index, row in cellproximities_concat_marshall.iterrows():

    cluster_A = row['unified_int'].split('--')[0]
    cluster_B = row['unified_int'].split('--')[1]

    if ((cluster_A, cluster_B) not in feasible_pairs_marshall)\
        |((cluster_A, cluster_B) not in feasible_pairs_marshall):
        feasible_pairs_marshall.append((cluster_A, cluster_B))

cellchat_marshall_wt = pd.concat(cellchat_marshall_wt_list)
cellchat_marshall_ki = pd.concat(cellchat_marshall_ki_list)
# Filter out insignificant communications
pval_threshold = 0.05
quantile = 0.25

cellchat_marshall_wt = cellchat_marshall_wt[(cellchat_marshall_wt['pval'] < pval_threshold)&(cellchat_marshall_wt['prob'] > cellchat_marshall_wt['prob'].quantile(quantile))]
cellchat_marshall_ki = cellchat_marshall_ki[(cellchat_marshall_ki['pval'] < pval_threshold)&(cellchat_marshall_ki['prob'] > cellchat_marshall_ki['prob'].quantile(quantile))]

# Construct the flow paths 
conditions = ['WT', 'KI']
cellchat_output_marshall = {'WT':cellchat_marshall_wt, 'KI':cellchat_marshall_ki}

# Construct the base network now
cfs.pp.construct_base_networks(adata_marshall,
                                cellchat_output_marshall,
                                condition_label='sample',
                                method='cellchat',
                                feasible_pairs=feasible_pairs_marshall,
                                celltype_sep_old='_',
                                celltype_sep_new='-',
                                node_sep='_',
                                base_network_label='base_networks')

# Construct the cell type ligand expression
cfs.pp.construct_celltype_ligand_expression(adata_marshall,
                                            celltype_label='cell_type',
                                            celltype_sep_old='_',
                                            celltype_sep_new='-',
                                            node_sep='_',
                                            expressions_label='X_celltype_ligand',
                                            base_network_label='base_networks')

# Learn the causal signaling networks now
cfs.tl.learn_causal_network(adata_marshall,
                            condition_label='sample',
                            control_label='WT',
                            causal_network_label='causal_networks',
                            celltype_ligands_label='X_celltype_ligand',
                            base_network_label='base_networks', 
                            n_cores=4,
                            n_bootstraps=100,
                            n_shuffles=1000,
                            alpha_ci=0.001,
                            alpha_inv=0.001)

# Validate causal network against base network
cfs.tl.validate_against_base_network(adata_marshall,
                            condition_label='sample',
                            control_label='WT',
                            causal_network_label='causal_networks',
                            celltype_ligands_label='X_celltype_ligand',
                            base_network_label='base_networks',
                            feasible_pairs=feasible_pairs_marshall,
                            celltype_sep_old='_',
                            celltype_sep_new='-')
