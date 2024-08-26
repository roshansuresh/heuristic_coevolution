# -*- coding: utf-8 -*-
"""
Pareto Fronts at different NFE - heuristic weights for a single run

@author: roshan94
"""
from Utils.dataHandler import DataHandler
import numpy as np
import os
import matplotlib.pyplot as plt

results_dir = 'C:\\SEAK Lab\\Coev results\\'
problem_dir = 'Assign\\' # Truss, Artery ,Assign, Partition, (C1_DTLZ1, Simple DTLZ1_2, UF1 are test problems)

int_pop_updated = False
if int_pop_updated:
    int_pop_dir = 'updated int pop\\'
else:
    int_pop_dir = 'same int pop\\'
    
int_weights = False # whether heuristic weights are integers or real values

obj_names = ['TrueObjective1','TrueObjective2']
if problem_dir == 'Truss\\' or problem_dir == 'Artery\\':
    heuristic_names = ['PartColl','NodalProp','Orient','Inters']
    #obj_names = ['Normalized Stiffness', 'Normalized Volume Fraction']
    # Set parameters for DataHandler.get_objectives() method
    objs_norm_num = [0, 0] 
    objs_norm_den = [1.8162e6, 1] # Youngs modulus used to normalize stiffness
    objs_max = [False, False] # both true objectives are stored so no need to multiply by -1
    if problem_dir == 'Artery\\':
        #obj_names = ['Normalized Stiffness', 'Normalized Deviation']
        # Set parameters for DataHandler.get_objectives() method
        objs_norm_num = [2e5, 0]
        objs_norm_den = [1e6, 1]
elif problem_dir == 'Assign\\':
    heuristic_names = ['Instrdc','Instrorb','Interinstr','Packeff','Spmass','Instrsyn','Instrcount']
    #obj_names = ['Normalized Science Score', 'Normalized Cost']
    # Set parameters for DataHandler.get_objectives() method
    objs_norm_num = [0, 0]
    objs_norm_den = [0.425, 2.5e4]
    objs_max = [True, False] # first objective (science) is to be maximized and second objective (cost) is to be minimized
else: # problem_dir == 'Partitioning Problem\\' (test problems not considered)
    heuristic_names = ['Instrdc','Instrorb','Interinstr','Packeff','Spmass','Instrsyn']
    # Set parameters for DataHandler.get_objectives() method
    objs_norm_num = [0, 0]
    objs_norm_den = [0.4, 7250]
    objs_max = [True, False] # first objective (science) is to be maximized and second objective (cost) is to be minimized
    
run_num = 0

heurs_incorporated = [True for i in range(len(heuristic_names))]
#heurs_incorporated = [False, False, True, True]

heurs_dir = ''
for i in range(len(heurs_incorporated)):
    if heurs_incorporated[i]:
        heurs_dir += heuristic_names[i]
        
heurs_dir += '\\'
current_filepath = results_dir + problem_dir + heurs_dir + int_pop_dir + 'run ' + str(run_num) + '\\'

# To determine internal population size, search the final coevolutionary population file and get the number of elements
internal_finalpop_filename = 'coevolutionary_algorithm_heuristic_weights_finalpop.csv'
internal_datahandler = DataHandler(current_filepath + internal_finalpop_filename)
internal_finalpop_columns = internal_datahandler.read(ignore_nans=False)
internal_pop_size = internal_datahandler.get_line_count()

# Read file contents and sort by NFE
coev_filename = 'coevolutionary_algorithm_heuristic_weights.csv'
coev_full_filename = current_filepath + coev_filename

data_handler = DataHandler(coev_full_filename)
file_columns = data_handler.read(ignore_nans=True)
sorted_file_columns = data_handler.sort_by_nfe()

coev_nfe_vals = sorted_file_columns.get('NFE')
pop_size = 0
for nfe in coev_nfe_vals:
    if nfe > 0:
        pop_size = nfe
        break

for i in range(int(pop_size)):
    coev_nfe_vals[i] = float(i)
max_coev_nfe = np.max(coev_nfe_vals)
heur_weights = data_handler.get_heur_weights()

# Finding best fitness values until current NFE
fitness = sorted_file_columns.get('Fitness Value 0')
best_fitness = np.zeros((len(fitness)))
best_fitness[0] = -1*fitness[0]
for i in range(1,len(coev_nfe_vals)):
    best_fitness[i] = np.max(np.multiply(fitness[:(i+1)], -1))

# Plotting fitness values 
fig_fit = plt.figure()
plt.plot(coev_nfe_vals, best_fitness, marker='*') # Negative of fitness values are stored since objective minimizer is used in the data generation script
plt.xlabel('NFE')
plt.ylabel('Coevolutionary Fitness')

if problem_dir == 'Equal Normal Stiffness Problem\\' or problem_dir == 'Artery Problem\\':
    num_feas = sorted_file_columns.get('Number of Feasible Solutions')
    fig_feas = plt.figure()
    plt.plot(coev_nfe_vals, num_feas, marker='*')
    plt.xlabel('NFE')
    plt.ylabel('Number of Internal Feasible Solutions')

# Search and read internal MOEA results for each set of heuristic weights
pareto_fronts = {}
recorded_coev_vals = []
for i in range(len(coev_nfe_vals)):
    nfe = coev_nfe_vals[i]
    heur_weight = heur_weights[i,:]
    
    # Construct heuristic weights string to check in filename
    check_str = ''
    for j in range(heur_weight.shape[0]):
        if int_weights:
            check_str += 'w' + str(j) + '-' + str(int(heur_weight[j])) + ';0' + '_'
        else:
            check_str += 'w' + str(j) + '-' + str(round(10**heur_weight[j], 5)).replace('.',';') + '_' # The values stored in the csv file are raw design decisions, must be converted to actual weights
    # Search in the current directory to find the results file for the given set of heuristic weights and in the current coev_nfe
    file_internal_MOEA = ''
    for root, dirs, files in os.walk(current_filepath):
        for filename in files:
            if (check_str in filename) and ('finalarchive' in filename): # check if filename contains the given heuristic weights
                # Next, check that the NFE value is close to `nfe' (due to limitations in the MOEA Framework implementation the saved NFE is the nearest previous population size multiple)
                str_nfe_start_ind = filename.find('nfe-')
                str_nfe_end_ind = filename.rfind('_')
                nfe_filename = int(filename[str_nfe_start_ind+4:str_nfe_end_ind])
                if np.absolute(nfe - nfe_filename) <= internal_pop_size:
                    file_internal_MOEA = filename
                    break
        break
    
    # Read the filename, extract and save Pareto objectives
    if file_internal_MOEA == '':
        continue # INVESTIGATE!! (possible reason: the same set of weights were evaluated either previously or subsequently)
    recorded_coev_vals.append(nfe)
    file_datahandler = DataHandler(current_filepath + file_internal_MOEA)
    internal_file_columns = file_datahandler.read(ignore_nans=False) # no need to sort, since this is the final archive
    pareto_objs = file_datahandler.get_objectives(obj_names=obj_names, objs_max=objs_max, objs_norm_den=objs_norm_den, objs_norm_num=objs_norm_num)
    pareto_fronts['Weights: ' + str(heur_weight) +', NFE:' + str(nfe)] = pareto_objs
    
# Plot Pareto Fronts for each coev NFE
n_figs = 5
n_inds = int(len(recorded_coev_vals)/n_figs)
pf_keys = list(pareto_fronts.keys())
if len(coev_nfe_vals) % n_figs > 0:
    n_figs += 1
current_nfe_idx = 0
for i in range(n_figs):
    fig = plt.figure()
    for j in range(n_inds):
        current_key = pf_keys[(n_inds*i)+j]
        str_nfe_start = current_key.find('NFE:')
        current_nfe = float(current_key[str_nfe_start+4:])
        pfs_nfe = pareto_fronts[current_key]
        
        # Assuming two objectives
        pfs_obj1 = [x[0] for x in pfs_nfe]
        pfs_obj2 = [x[1] for x in pfs_nfe]
        plt.scatter(pfs_obj1, pfs_obj2, marker='*', label='NFE = ' + str(current_nfe))
        plt.xlabel(obj_names[0])
        plt.ylabel(obj_names[1])
        plt.legend(loc='best')
    
    