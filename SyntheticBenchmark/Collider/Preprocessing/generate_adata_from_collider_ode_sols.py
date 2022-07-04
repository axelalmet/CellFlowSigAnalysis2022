import h5py
import numpy as np
import scanpy as sc
import pandas as pd

cwd = '../'

scenarios = ['collider']
num_simulations = 10

sol_columns = ['t', 'x', 'y', 'type', 'L1', 'R1', 'LR1', 'L2', 'R2', 'LR2', 'L3', 'R3', 'LR3']
var_names = ['CCL5', 'CCR1', 'GCG', 'GCGR', 'SST', 'SSTR1']
var_names_full = ['CCL5', 'CCR1', 'CCL5+CCR1', 'GCG', 'GCGR', 'GCG+GCGR', 'SST', 'SSTR1', 'SST+SSTR1']

for scenario in scenarios:   
    results_directory = cwd + 'output/' + scenario + '/'

    obs_data = []
    int_data = []

    obs_data_full = []
    int_data_full = []

    obs_results_directory = results_directory + 'observations/'
    int_results_directory = results_directory + 'interventions/'

    for seed in range(num_simulations):

        # Load the solutions
        obs_sol = pd.read_csv(obs_results_directory + scenario + "_causal_network_obs_sol_" + str(seed) + ".csv")
        int_sol = pd.read_csv(int_results_directory + scenario + "_causal_network_int_sol_" + str(seed) + ".csv")

        # Only keep the solutions at the end
        obs_sol = obs_sol[obs_sol['t'] == obs_sol['t'].max()]
        int_sol = int_sol[int_sol['t'] == int_sol['t'].max()]

        # Fix the columns (stupid formatting)
        obs_sol.columns = sol_columns
        int_sol.columns = sol_columns

        ### Construct the expressions

        # Ligand and total receptor
        obs_sol_expressions = np.zeros((obs_sol.shape[0], 6))
        obs_sol_expressions[:, 0] = obs_sol['L1']
        obs_sol_expressions[:, 1] = obs_sol['R1'] + obs_sol['LR1']
        obs_sol_expressions[:, 2] = obs_sol['L2']
        obs_sol_expressions[:, 3] = obs_sol['R2'] + obs_sol['LR2']
        obs_sol_expressions[:, 4] = obs_sol['L3']
        obs_sol_expressions[:, 5] = obs_sol['R3'] + obs_sol['LR3']

        int_sol_expressions = np.zeros((int_sol.shape[0], 6))
        int_sol_expressions[:, 0] = int_sol['L1']
        int_sol_expressions[:, 1] = int_sol['R1'] + int_sol['LR1']
        int_sol_expressions[:, 2] = int_sol['L2']
        int_sol_expressions[:, 3] = int_sol['R2'] + int_sol['LR2']
        int_sol_expressions[:, 4] = int_sol['L3']
        int_sol_expressions[:, 5] = int_sol['R3'] + int_sol['LR3']

        # Construct the annotated data frames
        adata_obs = sc.AnnData(X=obs_sol_expressions, obs=obs_sol[['x', 'y', 'type']], var=pd.DataFrame(index=var_names, data={'type':3*['ligand', 'receptor']}))
        adata_int = sc.AnnData(X=int_sol_expressions, obs=int_sol[['x', 'y', 'type']], var=pd.DataFrame(index=var_names, data={'type':3*['ligand', 'receptor']}))

        obs_data.append(adata_obs)
        int_data.append(adata_int)

        # Ligand, free receptor, and bound complex
        obs_sol_full_expressions = obs_sol[['L1', 'R1', 'LR1', 'L2', 'R2', 'LR2', 'L3', 'R3', 'LR3']].to_numpy()
        int_sol_full_expressions = int_sol[['L1', 'R1', 'LR1', 'L2', 'R2', 'LR2', 'L3', 'R3', 'LR3']].to_numpy()

        # Construct the annotated data frames
        adata_obs = sc.AnnData(X=obs_sol_expressions, obs=obs_sol[['x', 'y', 'type']], var=pd.DataFrame(index=var_names, data={'type':3*['ligand', 'receptor']}))
        adata_int = sc.AnnData(X=int_sol_expressions, obs=int_sol[['x', 'y', 'type']], var=pd.DataFrame(index=var_names, data={'type':3*['ligand', 'receptor']}))

        obs_data.append(adata_obs)
        int_data.append(adata_int)

        adata_obs_full = sc.AnnData(X=obs_sol_full_expressions, obs=obs_sol[['x', 'y', 'type']], var=pd.DataFrame(index=var_names_full, data={'type':3*['ligand', 'receptor', 'complex']}))
        adata_int_full = sc.AnnData(X=int_sol_full_expressions, obs=int_sol[['x', 'y', 'type']], var=pd.DataFrame(index=var_names_full, data={'type':3*['ligand', 'receptor', 'complex']}))

        obs_data_full.append(adata_obs_full)
        int_data_full.append(adata_int_full)

    # Ligand and total receptors
    adata_obs_joined = obs_data[0].concatenate(obs_data[1:], join='outer')
    adata_int_joined = int_data[0].concatenate(int_data[1:], join='outer')

    adata_obs_joined.obs['Condition'] = 'Observation'
    adata_obs_joined.obs_names_make_unique(join='_')

    adata_int_joined.obs['Condition'] = 'Intervention'
    adata_int_joined.obs_names_make_unique(join='_')

    adata_scenario = adata_obs_joined.concatenate(adata_int_joined, join='outer')
    adata_scenario.obs_names_make_unique(join='_')

    adata_scenario.obs['type'] = pd.Series(adata_scenario.obs['type'], dtype='category')

    # Rename the categories for CellChat
    adata_scenario.obs['type'].cat.rename_categories({0.0:'A', 1.0:'B', 2.0:'C', 3.0:'D', 4.0:'E'}, inplace=True)
    
    sc.pp.filter_cells(adata_scenario, min_counts=1)

    # Normalise and log-transform data and write it to file
    sc.pp.normalize_total(adata_scenario, target_sum=1e1)
    sc.pp.log1p(adata_scenario)

    adata_scenario.X[np.isnan(adata_scenario.X)] = 0.0

    adata_scenario.write(results_directory + 'adata_' + scenario + '.h5ad', compression='gzip')
