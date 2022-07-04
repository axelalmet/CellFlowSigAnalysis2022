### Script to simulate the evolution of a ligand-receptor model. As a first step, we're considering just one ligand, one receptor,
### and one complex on a network.
import networkx as nx
import matplotlib.pyplot as plt
import random
from scipy import integrate
import numpy as np
from itertools import count, product
from timeit import default_timer as timer
from numbalsoda import lsoda_sig, lsoda, address_as_void_pointer
from numba import njit, cfunc, carray, types
import numba as nb
import pickle

num_species = 3
num_celltypes = 5

# Some helper functions
def set_graph_attributes(G, model_params):

    # Get the parameters
    diffs = model_params['diffs']
    binds = model_params['binds']
    diss = model_params['diss']
    prods = model_params['prods']
    degrads = model_params['degrads']

    # Add the position attribute
    index = 0
    for node in G.nodes():
        i, j = node
        pos = np.array([10*i, 10*j])
        G.nodes[i, j]['pos'] = pos
        G.nodes[i, j]['index'] = index

        for k in range(num_species):
            G.nodes[i, j]['diff_' + str(k + 1)] = diffs[k] # Diffusion rate
            G.nodes[i, j]['bind_' + str(k + 1)]  = binds[k] # Binding rate
            G.nodes[i, j]['diss_' + str(k + 1)] = diss[k] # Dissociation rate
            G.nodes[i, j]['p_' + str(k + 1)]  = prods[k] # Auto-production term
            G.nodes[i, j]['d_' + str(k + 1)] = degrads[k] # Degradation rate

        # Let's also assign cell type while we're here
        if i < float(num_cols) / num_celltypes: # Cell type D is on the far left
            G.nodes[i, j]['type'] = 3
            
            # Set the production rates to zero for ligands 2 and 3 to be zero
            G.nodes[i, j]['p_1']  = 0
            G.nodes[i, j]['p_2']  = 0
            G.nodes[i, j]['p_3']  = 0

        elif ( (i >= float(num_cols) / num_celltypes)&(i < 2.0 * float(num_cols) / num_celltypes) ): # Cell type B is to the right of cell type D

            G.nodes[i, j]['type'] = 1

            # Set the production rates of Ligand 1 and 3 to be zero
            G.nodes[i, j]['p_1']  = 0
            G.nodes[i, j]['p_3']  = 0

        elif ( (i >= 2.0 * float(num_cols) / num_celltypes)&(i < 3.0 * float(num_cols) / num_celltypes) ): # Cell type A is in the middle

            G.nodes[i, j]['type'] = 0

            # Set the production rates of Ligand 1 and 2 to be zero
            G.nodes[i, j]['p_2']  = 0
            G.nodes[i, j]['p_3']  = 0

        elif ( (i >= 3.0 * float(num_cols) / num_celltypes)&(i < 4.0 * float(num_cols) / num_celltypes) ): # Cell type C is to the right of cell type A

            G.nodes[i, j]['type'] = 2

            # Set the production rates of Ligand 1 and 2 to be zero
            G.nodes[i, j]['p_1']  = 0
            G.nodes[i, j]['p_2']  = 0

        else: # Cell type E is on the far right

            G.nodes[i, j]['type'] = 4

            # No ligand production here
            G.nodes[i, j]['p_1']  = 0
            G.nodes[i, j]['p_2']  = 0
            G.nodes[i, j]['p_3']  = 0

        index += 1

def set_ode_parameters(G, initial_states, domain_height, domain_width, bws_ligand, bws_receptor, 
                        diff_rates, bind_rates, diss_rates, prod_rates, degrad_rates,
                        celltypes):

    num_nodes = len(G.nodes())

    for k in range(num_species):

        lig_x, lig_y, lig_r = np.random.rand(3)
        rec_x, rec_y, rec_r = np.random.rand(3)
        rec_x2, rec_y2, rec_r2 = np.random.rand(3)

        # Rescale the y position of the ligand and receptor
        lig_y = (domain_height) * lig_y
        rec_y = (domain_height) * rec_y

        rec_y2 = (domain_height) * rec_y2

        # Adjust the blob widths of the ligand and receptor
        lig_r = bws_ligand[0] + (bws_ligand[1] - bws_ligand[0])*lig_r
        rec_r = bws_receptor[0] + (bws_receptor[1] - bws_receptor[0])*rec_r

        rec_r2 = bws_receptor[0] + (bws_receptor[1] - bws_receptor[0])*rec_r2

        if k == 0:
            lig_x = (domain_width / num_celltypes) * (2.0 + lig_x)
            rec_x = (domain_width / num_celltypes) * (1.0 + rec_x)
            rec_x2 = (domain_width / num_celltypes) * (3.0 + rec_x2)

        elif k == 1:

            lig_x = (domain_width / num_celltypes) * (1.0 + lig_x)
            rec_x = (domain_width / num_celltypes) * rec_x

        else: # K = 2

            lig_x = (domain_width / num_celltypes) * (3.0 + lig_x)
            rec_x = (domain_width / num_celltypes) * (4.0 + rec_x)

        for node in G.nodes():
            i, j = node

            # Get the kinetic parameters for the ligand-receptor binding
            diff = G.nodes[i, j]['diff_' + str(k + 1)]
            bind_rate = G.nodes[i, j]['bind_' + str(k + 1)]
            diss_rate = G.nodes[i, j]['diss_' + str(k + 1)]
            prod_rate =  G.nodes[i, j]['p_' + str(k + 1)]
            degrad_rate =  G.nodes[i, j]['d_' + str(k + 1)]

            # Get the cell type and position
            celltype = G.nodes[i, j]['type']
            position = G.nodes[i, j]['pos']
            index = G.nodes[i, j]['index']
            
            diff_rates[k*num_nodes + index] = diff
            bind_rates[k*num_nodes + index] = bind_rate
            diss_rates[k*num_nodes + index] = diss_rate
            prod_rates[k*num_nodes + index] = prod_rate
            degrad_rates[k*num_nodes + index] = degrad_rate

            celltypes[index] = celltype

            # Set the ligand, receptor, and complex concentration for this node
            ligand = np.exp( -( (position[0] - lig_x)**2.0 + (position[1] - lig_y)**2.0 )/ (2.0 * lig_r)**2.0)
            receptor = np.exp( -( (position[0] - rec_x)**2.0 + (position[1] - rec_y)**2.0 )/ (2.0 * rec_r)**2.0)

            if k == 0:
                receptor += np.exp( -( (position[0] - rec_x2)**2.0 + (position[1] - rec_y2)**2.0 )/ (2.0 * rec_r2)**2.0)

            # Set the complex to be zero always
            complex = 0.0

            if k == 0: # For LR 1, cell type A is initialised with ligand 1, cell type B is initialised with receptor 2

                if celltype != 0: # Only cell type A has non-zero ligand
                    ligand = 0.0

                if ( (celltype != 1)&(celltype != 2) ): # Only cell type B has non-zero receptor
                    receptor = 0.0 # Set receptor to be zero

            elif k == 1: # For LR 2, cell type B is initialised with ligand 2, cell type C is initialised with receptor 2
                
                if celltype != 1: # Only cell type B has non-zero ligand
                    ligand = 0.0 # Set ligand to be zero

                if celltype != 3: # Only cell type C has non-zero receptor
                    receptor = 0.0 # Set receptor to be zero

            else: # Should be k = 2, i.r. LR 3, cell type C is initialised with ligand 3, cell type D is initialised with receptor 3

                if celltype != 2: # Only cell type C has non-zero ligand
                    ligand = 0.0 # Set ligand to be zero

                if celltype != 4: # Only cell type D has non-zero receptor
                    receptor = 0.0 # Set receptor to be zero

            # We should be able to set the initial state now

            initial_states[3*k*num_nodes + 3*index] = ligand
            initial_states[3*k*num_nodes + 3*index + 1] = receptor
            initial_states[3*k*num_nodes + 3*index + 2] = complex

            G.nodes[i, j]['L_' + str(k + 1)] = ligand
            G.nodes[i, j]['R_' + str(k + 1)] = receptor
            G.nodes[i, j]['LR_' + str(k + 1)] = complex

def set_ode_sol_at_time(G, time, all_timepoints, sol):

    timeIndex = np.where(all_timepoints == time)[0][0]  
    node_pairs = list(G.nodes())

    for k in range(num_species):
        for node in node_pairs:
            i, j = node
            index = G.nodes[i, j]['index']

            G.nodes[i, j]['L_' + str(k + 1)] = sol[timeIndex, 3*k*num_nodes + 3*index]
            G.nodes[i, j]['R_' + str(k + 1)] = sol[timeIndex, 3*k*num_nodes + 3*index + 1]
            G.nodes[i, j]['LR_' + str(k + 1)] = sol[timeIndex, 3*k*num_nodes + 3*index + 2]

def plot_lr_species(G, fig_width, fig_height, v_min, v_max):
   
    plot_count = 1
    node_coordinates = dict(G.nodes(data='pos')) # Store as node coordinates for later

    # Initialise figure
    plt.figure(figsize=(fig_width, fig_height))

    for k in range(num_species):
        plt.subplot(num_species, 3, plot_count)
        choice = 'L_' + str(k + 1)
        node_values = set(nx.get_node_attributes(G, choice).values())
        mapping = dict(zip(sorted(node_values), count()))
        nodes = G.nodes()
        colours = [mapping[G.nodes[node][choice]] for node in nodes]
        colours_norm = plt.Normalize(min(colours), max(colours))
        ec = nx.draw_networkx_edges(G, node_coordinates, alpha = 0.2)
        nc = nx.draw_networkx_nodes(G, node_coordinates, nodelist=nodes, node_color = colours, node_size=10, cmap = plt.cm.Spectral_r, vmin=v_min, vmax=v_max)
        plt.colorbar(nc)
        plt.axis('off')
        plot_count += 1 # Update plot count

        plt.subplot(num_species, 3, plot_count)
        choice = 'R_' + str(k + 1)
        node_values = set(nx.get_node_attributes(G, choice).values())
        mapping = dict(zip(sorted(node_values), count()))
        nodes = G.nodes()
        colours = [mapping[G.nodes[node][choice]] for node in nodes]
        ec = nx.draw_networkx_edges(G, node_coordinates, alpha = 0.2)
        nc = nx.draw_networkx_nodes(G, node_coordinates, nodelist=nodes, node_color = colours, node_size=10, cmap = plt.cm.Spectral_r, vmin=v_min, vmax=v_max)
        plt.colorbar(nc)
        plt.axis('off')
        plot_count += 1 # Update plot count

        plt.subplot(num_species, 3, plot_count)
        choice = 'LR_' + str(k + 1)
        node_values = set(nx.get_node_attributes(G, choice).values())
        mapping = dict(zip(sorted(node_values), count()))
        nodes = G.nodes()
        colours = [mapping[G.nodes[node][choice]] for node in nodes]
        ec = nx.draw_networkx_edges(G, node_coordinates, alpha = 0.2)
        nc = nx.draw_networkx_nodes(G, node_coordinates, nodelist=nodes, node_color = colours, node_size=10, cmap = plt.cm.Spectral_r, vmin=v_min, vmax=v_max)
        # nc = nx.draw_networkx_nodes(G, node_coordinates, nodelist = nodes, node_size=100)
        plt.colorbar(nc)
        plt.axis('off')
        plot_count += 1 # Update plot count

    plt.show() 

def construct_ode_sol(G, times, ode_sol):

    num_times = len(times)
    
    node_pairs = list(G.nodes())
    num_nodes = len(node_pairs)

    # Number of columns comes from: time, x, y, L1, R1, LR1, ..., Ln, Rn, LRn
    output_sol = np.zeros((num_times*num_nodes, 4 + 3*num_species))

    # Set positions
    for t in range(num_times):

        time = times[t]

        for k in range(num_species):

            for node in G.nodes():

                i, j = node
                index = G.nodes[i, j]['index']
                celltype = G.nodes[i, j]['type']
                position = G.nodes[i, j]['pos']

                output_sol[t*num_nodes + index, 0] = time
                output_sol[t*num_nodes + index, 1] = position[0].copy()
                output_sol[t*num_nodes + index, 2] = position[1].copy()
                output_sol[t*num_nodes + index, 3] = celltype
                output_sol[t*num_nodes + index, 4 + 3*k] = ode_sol[t, 3*k*num_nodes + 3*index]
                output_sol[t*num_nodes + index, 4 + 3*k + 1] = ode_sol[t, 3*k*num_nodes + 3*index + 1]
                output_sol[t*num_nodes + index, 4 + 3*k + 2] = ode_sol[t, 3*k*num_nodes + 3*index + 2]
    
    return output_sol

### Set-up for analysis
simulation_path = '../output/fork/'
observational_path = simulation_path + 'observations'
interventional_path = simulation_path + 'interventions'

# Define the causality links
causal_prods = np.zeros((num_species, num_species), dtype=np.float64)
causal_prods[0, 1] = 0.1 # LR complex 1 triggers L2 
causal_prods[0, 2] = 0.1 # LR complex 2 triggers L3

celltype_causality = np.zeros((num_species, num_celltypes), dtype=np.float64)
celltype_causality[0, 1] = 1.0 # LR complex 1 triggers cell type 2
celltype_causality[0, 2] = 1.0 # LR complex 2 triggers cell type 3

# Define the parameters
diffs = np.array([1.0, 1.0, 1.0]) # Diffusivity of ligands
binds = np.array([0.2, 0.2, 0.2]) # Binding rate of ligand-receptor complexes
diss = np.array([0.0, 0.0, 0.0]) # Dissociation rates of ligand-receptor complexes
degrads = np.array([0.01, 0.01, 0.01]) # Degradation rate (2 for each pair)
prods = np.array([0.1, 0.1, 0.1])

# Generate a 2D point
num_rows = 10
num_cols = 50
num_nodes = num_rows*num_cols # Total number of nodes
domain_width = 10.0*num_cols
domain_height = 10.0*num_rows

# Generate 2D grid with 4 nearest neighbours
G = nx.grid_2d_graph(num_cols, num_rows)
D = -1.0 * nx.laplacian_matrix(G).toarray().astype('float64') # Define the Laplacian matrix.
A = nx.adjacency_matrix(G) # Adjacency

model_params = {'diffs':diffs, 'binds':binds, 'diss':diss, 'degrads':degrads, 'prods':prods,
                 'causal_prods':causal_prods, 'celltype_causality':celltype_causality,
                 'num_rows':num_rows, 'num_cols':num_cols, 'A':A, 'D':D}

# Set all the graph attributes
set_graph_attributes(G, model_params)

# Define the node pairs
node_pairs = list(G.nodes())

# Collect the different kinetic parameters
diff_rates = np.zeros(num_species*num_nodes, dtype=np.float64)
bind_rates = np.zeros(num_species*num_nodes, dtype=np.float64)
diss_rates = np.zeros(num_species*num_nodes, dtype=np.float64)
prod_rates = np.zeros(num_species*num_nodes, dtype=np.float64)
degrad_rates = np.zeros(num_species*num_nodes, dtype=np.float64)
celltypes = np.zeros(num_nodes, dtype=int)

for k in range(num_species):
    for node in node_pairs:
        i, j = node
        index = G.nodes[i, j]['index']
        cell_type = G.nodes[i, j]['type']
        celltypes[index] = cell_type

        diff_rates[k*num_nodes + index] = G.nodes[i, j]['diff_' + str(k + 1)] # Diffusion rate
        bind_rates[k*num_nodes + index] = G.nodes[i, j]['bind_' + str(k + 1)] # Binding rate
        diss_rates[k*num_nodes + index] = G.nodes[i, j]['diss_' + str(k + 1)] # Dissociation rate
        prod_rates[k*num_nodes + index] = G.nodes[i, j]['p_' + str(k + 1)] # Auto-production term
        degrad_rates[k*num_nodes + index] = G.nodes[i, j]['d_' + str(k + 1)] # Degradation rate

# Generate the initial state
blobwidths_ligand = [10.0, 20.0]
blobwidths_receptor = [20.0, 40.0]

### Define the ODE function
def fork_rhs(t, Y, dY, D, causal_prods, celltype_causality):
    
    for k in range(num_species):

        for index in range(num_nodes):

            # Get the kinetic parameters for the ligand-receptor binding
            diff = diff_rates[k*num_nodes + index]
            bind_rate = bind_rates[k*num_nodes + index]
            diss_rate = diss_rates[k*num_nodes + index]
            prod_rate = prod_rates[k*num_nodes + index]
            degrad_rate = degrad_rates[k*num_nodes + index]

            # Get the cell type
            celltype = celltypes[index]

            ligand = Y[3*(k*num_nodes + index)]
            receptor = Y[3*(k*num_nodes + index) + 1]
            complex = Y[3*(k*num_nodes + index) + 2]
            ligand_flux = sum(np.array([D[index, n]*Y[3*(k*num_nodes + n)] for n in range(num_nodes)]))

            causal_prod_term = sum(np.array([causal_prods[m, k]*celltype_causality[m, celltype]*Y[3*m*num_nodes + 3*index + 2] for m in range(num_species)]))

            # Dynamics of ligand: diffusion, binding, dissociation, production, and degradation
            dY[3*(k*num_nodes + index)] = diff * ligand_flux\
                                            - bind_rate * ligand * receptor\
                                            + diss_rate * complex\
                                            + prod_rate\
                                            - degrad_rate * ligand\
                                            + causal_prod_term
            dY[3*(k*num_nodes + index) + 1] = diss_rate * complex - bind_rate * ligand * receptor
            dY[3*(k*num_nodes + index) + 2] = bind_rate * ligand * receptor - diss_rate * complex

# Define the argument types of the different parameters. We need to specify the array types and the shape parametesr too
args_dtype = types.Record.make_c_struct([
                    ('D_p', types.int64),
                    ('causal_prods_p', types.int64),
                    ('celltype_causality_p', types.int64),
                    ('D_shape_0', types.int64),
                    ('D_shape_1', types.int64),
                    ('causal_prods_shape_0', types.int64),
                    ('causal_prods_shape_1', types.int64),
                    ('celltype_causality_shape_0', types.int64),
                    ('celltype_causality_shape_1', types.int64)])

# this function will create the numba function to pass to lsoda.
def create_jit_rhs(rhs, args_dtype):
    jitted_rhs = njit(rhs)
    @nb.cfunc(types.void(types.double,
             types.CPointer(types.double),
             types.CPointer(types.double),
             types.CPointer(args_dtype)))
    def wrapped(t, u, du, user_data_p):
        # unpack p and arr from user_data_p
        user_data = nb.carray(user_data_p, 1)
        D_p = nb.carray(address_as_void_pointer(user_data[0].D_p),(user_data[0].D_shape_0, user_data[0].D_shape_1), dtype=np.float64)
        caus_prods = nb.carray(address_as_void_pointer(user_data[0].causal_prods_p),(user_data[0].causal_prods_shape_0, user_data[0].causal_prods_shape_1), dtype=np.float64)
        celltype_caus = nb.carray(address_as_void_pointer(user_data[0].celltype_causality_p),(user_data[0].celltype_causality_shape_0, user_data[0].celltype_causality_shape_1), dtype=np.float64)
        
        # then we call the jitted rhs function, passing in data
        jitted_rhs(t, u, du, D_p, caus_prods, celltype_caus) 
    return wrapped

# JIT the ODE RHS
rhs_cfunc = create_jit_rhs(fork_rhs, args_dtype)

t_eval = np.linspace(0, 10, 2)
intervention_scale = 1.0

num_realisations = 10

# np.random.seed(1) # Set a random seed for reproducibility
for seed in range(num_realisations):
    
    initial_states = np.zeros(3*num_species*num_nodes)
    
    set_ode_parameters(G, initial_states, domain_height, domain_width, blobwidths_ligand, blobwidths_receptor, 
                            diff_rates, bind_rates, diss_rates, prod_rates, degrad_rates,
                            celltypes)

    args_obs = np.array((D.ctypes.data,
                causal_prods.ctypes.data,
                celltype_causality.ctypes.data,
                D.shape[0], D.shape[1],
                causal_prods.shape[0], causal_prods.shape[1],
                celltype_causality.shape[0], celltype_causality.shape[1]),dtype=args_dtype)

    funcptr = rhs_cfunc.address
    obs_sol, success = lsoda(funcptr, initial_states, t_eval, data=args_obs)

    # Save the model parameters
    obs_params = model_params.copy() 
    f = open('%s/meta_data_' % observational_path + str(seed) + '.pkl',"wb")
    pickle.dump(obs_params, f)
    f.close()

    # Save the ODE sol
    obs_sol_output = construct_ode_sol(G, t_eval, obs_sol)
    np.savetxt(observational_path + "/fork_causal_network_obs_sol_" + str(seed) + ".csv", obs_sol_output, delimiter=',', header="t,x,y,type,L1,R1,LR1,L2,R2,LR2,L3,R3,LR3", comments="")

    # Set new initial state
    set_ode_parameters(G, initial_states, domain_height, domain_width, blobwidths_ligand, blobwidths_receptor, 
                        diff_rates, bind_rates, diss_rates, prod_rates, degrad_rates,
                        celltypes)

    # Simulate a soft intervention by perturbing the causal production rates
    intervened_causal_prods = causal_prods.copy()

    # intervened_causal_prods[0, 1] += np.random.normal(loc=0.0, scale=intervention_scale)
    # intervened_causal_prods[0, 2] += np.random.normal(loc=0.0, scale=intervention_scale)
    intervened_causal_prods[0, 1] += np.random.lognormal(mean=0.0, sigma=intervention_scale)
    intervened_causal_prods[0, 2] += np.random.lognormal(mean=0.0, sigma=intervention_scale)

    # Define the new arguments for the interventional parameters
    args_int = np.array((D.ctypes.data,
                    intervened_causal_prods.ctypes.data,
                    celltype_causality.ctypes.data,
                    D.shape[0], D.shape[1],
                    intervened_causal_prods.shape[0], intervened_causal_prods.shape[1],
                    celltype_causality.shape[0], celltype_causality.shape[1]),dtype=args_dtype)

    funcptr = rhs_cfunc.address
    int_sol, success = lsoda(funcptr, initial_states, t_eval, data=args_int)

    # Save the model parameters
    int_params = model_params.copy() 
    int_params['causal_prods'] = intervened_causal_prods
    f = open('%s/meta_data_' % interventional_path + str(seed) + '.pkl',"wb")
    pickle.dump(int_params, f)
    f.close()

    int_sol_output = construct_ode_sol(G, t_eval, int_sol)
    np.savetxt(interventional_path + "/fork_causal_network_int_sol_" + str(seed) + ".csv", int_sol_output, delimiter=',', header="t,x,y,type,L1,R1,LR1,L2,R2,LR2,L3,R3,LR3", comments="")
