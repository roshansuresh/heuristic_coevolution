# -*- coding: utf-8 -*-
"""
Plot Fitness Stats across all runs for a single case

@author: roshan94
"""
from Utils.dataHandler import DataHandler
#import statistics
import numpy as np
import matplotlib.pyplot as plt

results_dir = 'C:\\SEAK Lab\\Coev results\\'
problem_dir = 'Partition\\' # Equal Stiffness, Artery, Assign, Partition, (C1_DTLZ1, Simple DTLZ1_2, UF1 are test problems)
        
weight_of_weights = False # whether an additional weight of weights design decision is used
period_zero_inj = False # whether the zero solution is injected into the population at each 

obj_names = ['TrueObjective1','TrueObjective2']
if problem_dir == 'Truss\\' or problem_dir == 'Artery\\':
    heuristic_names = ['P','N','O','I']
    #heuristic_names = ['PartColl','NodalProp','Orient','Inters']
    #obj_names = ['Normalized Stiffness', 'Normalized Volume Fraction']
    #if problem_dir == 'Artery\\':
        #obj_names = ['Normalized Stiffness', 'Normalized Deviation']
elif problem_dir == 'Assign\\':
    heuristic_names = ['D','O','I','P','M','S','C']
    #heuristic_names = ['Instrdc','Instrorb','Interinstr','Packeff','Spmass','Instrsyn','Instrcount']
    #obj_names = ['Normalized Science Score', 'Normalized Cost']
else: # problem_dir == 'Partitioning Problem\\' (test problems not considered)
    heuristic_names = ['D','O','I','P','M','S']
    #heuristic_names = ['Instrdc','Instrorb','Interinstr','Packeff','Spmass','Instrsyn']

n_runs = 30
include_init_pop = True
include_labels = True
coev_pop_size = 5 # Taken from java code 

case_booleans = [True for i in range(len(heuristic_names))] # All heuristics enforced
#case_booleans = [False, False, True, True]

fitness_case = {}

coev_dir = ''
if weight_of_weights:
    coev_dir = 'WoW\\'
elif period_zero_inj: 
    coev_dir = 'PZI\\'    
else:
    coev_dir = 'Base\\'    

heurs_dir = ''
for i in range(len(case_booleans)):
    if case_booleans[i]:
        heurs_dir = 'Coevolutionary\\'
        break

for i in range(n_runs):
    current_filepath = results_dir + problem_dir + heurs_dir + coev_dir + 'run ' + str(i) + '\\'
    
    # Read file contents and sort by NFE
    coev_filename = 'coevolutionary_algorithm_heuristic_weights.csv'
    coev_full_filename = current_filepath + coev_filename

    data_handler = DataHandler(coev_full_filename)
    file_columns = data_handler.read(ignore_nans=False)
    sorted_file_columns = data_handler.sort_by_nfe()

    coev_nfe_vals = sorted_file_columns.get('NFE')
    max_coev_nfe = np.max(coev_nfe_vals)
    heur_weights = data_handler.get_heur_weights()

    # Storing fitness values 
    fitness = sorted_file_columns.get('Fitness Value 0')
    fitness_case['run'+str(i)] = np.multiply(fitness,-1) # since Eps MOEA uses a minimizer, -fitness was used as objective for the weights solution 
    
# Compute fitness stats at each NFE
run_keys = fitness_case.keys()

n_nfes = len(coev_nfe_vals)

# Plot fitness stats
n_figs = 6
img_counter = 0
for i in range(n_figs):
    #plt.subplot(2,3,i+1)
    plt.figure()
    for j in range(int(n_runs/n_figs)):
        if include_init_pop:
            if include_labels:
                plt.plot(coev_nfe_vals, fitness_case['run'+str(img_counter)], label='run'+str(img_counter+1))
            else:
                plt.plot(coev_nfe_vals, fitness_case['run'+str(img_counter)])
        else:
            if include_labels:
                plt.plot(coev_nfe_vals[coev_pop_size:], fitness_case['run'+str(img_counter)][coev_pop_size:], label='run'+str(img_counter+1))
            else:
                plt.plot(coev_nfe_vals[coev_pop_size:], fitness_case['run'+str(img_counter)][coev_pop_size:])
        img_counter += 1
    plt.xlabel(r'Number of Function Evaluations',fontsize=12)
    plt.ylabel(r'Coevolutionary Fitness',fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    if include_labels:
        plt.legend(loc='best')
    
