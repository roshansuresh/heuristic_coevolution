# -*- coding: utf-8 -*-
"""
Pareto Fronts at different NFE - heuristic weights for a single run

@author: roshan94
"""
from Utils.dataHandler import DataHandler
from Utils.mOOCaseStatistics import MOOCaseStatistics
from Utils.mOORunStatistics import MOORunStatistics
import numpy as np
import os
import matplotlib.pyplot as plt

results_dir = 'C:\\SEAK Lab\\Coev results\\'
problem_dir = 'Artery\\' # Truss, Artery, Assign, Partition, (C1_DTLZ1, Simple DTLZ1_2, UF1 are test problems)

int_pop_updated = True # whether updated internal population is used in successive heuristic weight evaluations or not
ga_alg = True # Genetic Algorithm (GA) results if True, Differential Evolution (DE) results if False
constr_viol_fitness = True # Constraint Violation based weights fitness if True, feasible HV difference based weights fitness if False
int_weights = False # whether heuristic weights are integers or real values
weight_of_weights = False # whether an additional weight of weights design decision is used
period_zero_inj = False # whether the zero solution is injected into the population at each 
constr_names = ''  

obj_names = ['TrueObjective1','TrueObjective2']
if problem_dir == 'Truss\\' or problem_dir == 'Artery\\':
    heuristic_names = ['P','N','O','I']
    #heuristic_names = ['PartColl','NodalProp','Orient','Inters']
    #obj_names = ['Normalized Stiffness', 'Normalized Volume Fraction']
    # Set parameters for DataHandler.get_objectives() method
    objs_norm_num = [0, 0] 
    objs_norm_den = [1.8162e6, 1] # Youngs modulus used to normalize stiffness
    objs_max = [False, False] 
    # To be set to true if negative of any objective is to be used to compute HV, 
    # first objective (stiffness) is to be maximized and second objective (volume fraction/deviation) is to be minimized, however -normalized stiffness is stored in csv so -1 multiplication is not required
    true_obj_names = [r'$C_{22}$',r'$v_f$'] 
    if problem_dir == 'Artery\\':
        #obj_names = ['Normalized Stiffness', 'Normalized Deviation']
        # Set parameters for DataHandler.get_objectives() method
        objs_norm_num = [2e5, 0]
        objs_norm_den = [1e6, 1]
        true_obj_names = [r'$\frac{C_{11}}{v_f}$',r'deviation']
elif problem_dir == 'Assign\\':
    #heuristic_names = ['Instrdc','Instrorb','Interinstr','Packeff','Spmass','Instrsyn','Instrcount']
    heuristic_names = ['D','O','I','P','M','S','C']
    #obj_names = ['Normalized Science Score', 'Normalized Cost']
    # Set parameters for DataHandler.get_objectives() method
    objs_norm_num = [0, 0]
    objs_norm_den = [0.425, 2.5e4]
    objs_max = [False, False] 
    # To be set to true if negative of any objective is to be used to compute HV, 
    # first objective (science) is to be maximized and second objective (cost) is to be minimized, however -normalized science is stored in csv so -1 multiplication is not required
    true_obj_names = [r'Science',r'Cost']
else: # problem_dir == 'Partitioning Problem\\' (test problems not considered)
    #heuristic_names = ['Instrdc','Instrorb','Interinstr','Packeff','Spmass','Instrsyn']
    heuristic_names = ['D','O','I','P','M','S']
    # Set parameters for DataHandler.get_objectives() method
    objs_norm_num = [0, 0]
    objs_norm_den = [0.4, 7250]
    objs_max = [False, False] 
    # To be set to true if negative of any objective is to be used to compute HV, 
    # first objective (science) is to be maximized and second objective (cost) is to be minimized, however -normalized science is stored in csv so -1 multiplication is not required
    true_obj_names = [r'Science',r'Cost']
    
n_runs = 20

case_bools = {} # Last boolean signifies whether coevolutionary algorithm is used (as opposed to AOS)
case_bools['Eps. MOEA'] = [False for i in range(len(heuristic_names) + 1)] # No heuristics enforced
aos_case_bools = [True for i in range(len(heuristic_names))]
aos_case_bools.append(False)
case_bools['AOS'] =  aos_case_bools # All heuristics enforced through AOS
case_bools['Coev. Penalty'] = [True for i in range(len(heuristic_names) + 1)] # All heuristics enforced through Coevolutionary Interior Penalization

plot_colors = {} # black, red, blue for Eps. MOEA, AOS & Coev. Penalty respectively
plot_colors['Eps. MOEA'] = 'k'
plot_colors['AOS'] = 'r'
plot_colors['Coev. Penalty'] = 'b'

plot_markers = {}
plot_markers['Eps. MOEA'] = '*'
plot_markers['AOS'] = 'o'
plot_markers['Coev. Penalty'] = '+'

## Sample NFEs
n_samples = 2

######## POTENTIAL: Script to ascertain pop size from csv file

pop_size = 90
max_nfe = 13500
if problem_dir == 'Artery\\' or problem_dir == 'Truss\\':
    max_nfe = 13500
    
#nfe_samples = np.linspace(pop_size, max_nfe, n_samples)
nfe_samples = [pop_size, 270]
#nfe_samples = [pop_size, 300]

## Create dictionary of all solution objectives (and constraints if applicable) for all runs for each case (mainly collate internal runs for each coevolutionary run)
case_objs = {}
case_constrs = {}
case_nfes = {}
int_pop_found = False
print('Reading and storing solutions for all cases')

for case_key in list(case_bools.keys()):
    run_objs = {}
    run_constrs = {}
    run_nfes = {}
    heurs_incorporated = case_bools[case_key]
    print('Reading Case: ' + case_key)
    
    if int_pop_updated:
        int_pop_dir = 'updated int pop\\'
    else:
        int_pop_dir = 'same int pop\\'
        
    alg_dir = 'GA\\' 
    if not ga_alg: 
        alg_dir = 'DE\\'
        
    fitness_dir = ''
    if problem_dir == 'Truss\\' or problem_dir == 'Artery\\':
        fitness_dir = 'constraint violation fitness\\'
        if not constr_viol_fitness:
            fitness_dir = 'feasible hypervolume fitness\\'
            
    wow_dir = ''
    if weight_of_weights:
        wow_dir = 'WoW - alpha 10e-x\\'
        
    pzi_dir = ''
    if period_zero_inj:
        pzi_dir = 'PZI - alpha 10e-x\\'
    
    heurs_dir = ''
    for i in range(len(heurs_incorporated) - 1):
        if heurs_incorporated[i]:
            heurs_dir += heuristic_names[i]
            
    if heurs_dir == '':
        heurs_dir = 'Eps MOEA'
        int_pop_dir = ''
    elif (not heurs_incorporated[-1]):
        heurs_dir = 'AOS'
        int_pop_dir = ''
            
    heurs_dir += '\\'
    for i in range(n_runs):
        nfes_current_run = []
        current_run_objs = []
        current_run_constrs = []
        print('Reading Run: ' + str(i))
        if (heurs_dir == 'Eps MOEA\\') or (heurs_dir == 'AOS\\'): # Epsilon MOEA or AOS - just read csv file for corresponding run
            current_filepath = results_dir + problem_dir + heurs_dir + int_pop_dir
            
            if heurs_dir == 'Eps MOEA\\':
                run_filename = 'EpsilonMOEA_' + str(i) + '_allSolutions.csv'
            else:
                run_filename = 'AOSMOEA_' + str(i) + '_allSolutions.csv'
            
            runfile_datahandler = DataHandler(current_filepath + run_filename)
            runfile_cols = runfile_datahandler.read(ignore_nans=False)
            sorted_runfile_cols = runfile_datahandler.sort_by_nfe()
            
            nfes_current_run = sorted_runfile_cols.get('NFE')
            current_run_objs = runfile_datahandler.get_objectives(obj_names=obj_names, objs_max=objs_max, objs_norm_den=objs_norm_den, objs_norm_num=objs_norm_num)
            if not constr_names == '':
                current_run_constrs = runfile_datahandler.get_parameters(parameter_names=constr_names)
        
        else: # Coevolutionary MOEA for heuristic enforcement - read csv files associated with internal MOEA runs and store collated objectives and constraints
            current_filepath = results_dir + problem_dir + heurs_dir + int_pop_dir + alg_dir + fitness_dir + wow_dir + pzi_dir + 'run ' + str(i) + '\\'
            
            # Read file contents and sort by NFE
            coev_filename = 'coevolutionary_algorithm_heuristic_weights.csv'
            coev_full_filename = current_filepath + coev_filename
            
            data_handler = DataHandler(coev_full_filename)
            file_columns = data_handler.read(ignore_nans=True)
            sorted_file_columns = data_handler.sort_by_nfe()
            
            coev_nfe_vals = sorted_file_columns.get('NFE')
            heur_weights = data_handler.get_heur_weights()
            
            #nfes_runfile = []
            #objs_runfile = []
            #constrs_runfile = []
            for n in range(len(coev_nfe_vals)):
                coev_nfe = coev_nfe_vals[n]
                heur_weight = heur_weights[n,:]
                
                start_idx = 0
                if weight_of_weights:
                    start_idx = 1
                
                check_str = ''
                weight_mult = 1
                if weight_of_weights:
                    weight_mult = heur_weight[0]
                for j in range(len(heuristic_names)):
                    if int_weights:
                        check_str += 'w' + str(j) + '-' + str(int(heur_weight[j+start_idx])*weight_mult) + ';0' + '_'
                    else:
                        weight_val = round((10**heur_weight[j+start_idx])*weight_mult, 5) # The values stored in the csv file are raw design decisions, must be converted to actual weights
                        if weight_val < 1e-3 and weight_val > 0:
                            check_str += 'w' + str(j) + '-' + np.format_float_scientific(weight_val, exp_digits=1, precision=1, min_digits=1).upper().replace('.',';') + '_' 
                        else:
                            check_str += 'w' + str(j) + '-' + str(weight_val).replace('.',';') + '_' 
                
                if weight_of_weights:
                    if weight_mult < 1e-3 and weight_mult > 0:
                        check_str += 'ww-' + np.format_float_scientific(weight_mult, exp_digits=1, precision=1, min_digits=1).upper().replace('.',';') + '_'
                    else:
                        check_str += 'ww-' + str(round(weight_mult, 5)).replace('.',';') + '_'
                
                # Search in the current directory to find the results file for the given set of heuristic weights and in the current coev_nfe
                file_internal_MOEA = ''
                for root, dirs, files in os.walk(current_filepath):
                    for filename in files:
                        if (check_str in filename) and ('allSolutions' in filename): # check if filename contains the given heuristic weights
                            # Next, check that the NFE value is close to `nfe' (due to limitations in the MOEA Framework implementation the saved NFE is the nearest previous population size multiple)
                            str_nfe_start_ind = filename.find('nfe-')
                            str_nfe_end_ind = filename.rfind('_')
                            nfe_filename = int(filename[str_nfe_start_ind+4:str_nfe_end_ind])
                            if coev_nfe == nfe_filename:
                                file_internal_MOEA = filename
                                break
                    break
                
                # Read the filename, extract and save objectives and constraints
                file_datahandler = DataHandler(current_filepath + file_internal_MOEA)
                
                if not int_pop_found:
                    internal_file_allcols = file_datahandler.read(ignore_nans=True)
                    internal_pop_size = file_datahandler.get_line_count()
                    file_datahandler.reset()
                    int_pop_found = True
                
                internal_file_cols = file_datahandler.read(ignore_nans=False)
                internal_file_sorted_cols = file_datahandler.sort_by_nfe()
                internal_file_nfes = internal_file_sorted_cols.get('NFE')
                internal_objs = file_datahandler.get_objectives(obj_names=obj_names, objs_max=objs_max, objs_norm_den=objs_norm_den, objs_norm_num=objs_norm_num)
                if not constr_names == '':
                    internal_constr = file_datahandler.get_parameters(parameter_names=constr_names)                
                
                for j in range(len(internal_file_nfes)):
                    nfes_current_run.append(coev_nfe*internal_pop_size + internal_file_nfes[j])
                    current_run_objs.append(internal_objs[j,:])
                    if not constr_names == '':
                        current_run_constrs.append(internal_constr[j,:])
                    
        run_nfes['run'+str(i)] = nfes_current_run
        run_objs['run'+str(i)] = np.stack(current_run_objs, axis=0)
        if not constr_names == '':
            run_constrs['run'+str(i)] = np.stack(current_run_constrs, axis=0)
        
    case_nfes[case_key] = run_nfes
    case_objs[case_key] = run_objs
    if not constr_names == '':
        case_constrs[case_key] = run_constrs

## Create dictionary of Pareto Fronts at sample nfe for all runs for each case
print('Computing and storing Pareto Fronts')
dummy_caseStats = MOOCaseStatistics(hv_allcases={}, nfe_array=nfe_samples, case_names=list(case_bools.keys()))
case_pfs = {}
for case_key in list(case_objs.keys()):
    print('Computing for Case: ' + case_key)
    objs_case = case_objs[case_key]
    nfes_case = case_nfes[case_key]
    if not constr_names == '':
        constrs_case = case_constrs[case_key]
    #run_pfs = {}
    nfe_sample_pfs = {}
    for i in range(len(nfe_samples)):
        nfe_sample = nfe_samples[i]
        print('Computing for NFE: ' + str(nfe_sample))
        objs_runs_nfe = []
        constr_runs_nfe = []
        
        for run_key in list(objs_case.keys()):
            objs_run = objs_case[run_key]
            if not constr_names == '':
                constrs_run = constrs_case[run_key]
            
            nfes_run = nfes_case[run_key]
            
            pfs_run = {}
            objs_fullsat = []
            
            closest_nfe_idx = dummy_caseStats.find_closest_index(val=nfe_sample, search_list=nfes_run)
            objs_nfe_sample = objs_run[:closest_nfe_idx,:]
            #objs_nfe_unique, nfe_unique_idx = np.unique(objs_nfe_sample, return_index=True, axis=0)
            
            if len(objs_runs_nfe) == 0:
                #objs_runs_nfe = objs_nfe_unique
                objs_runs_nfe = objs_nfe_sample
            else:
                #objs_runs_nfe = np.append(objs_runs_nfe, objs_nfe_unique, axis=0)
                objs_runs_nfe = np.append(objs_runs_nfe, objs_nfe_sample, axis=0)
            
            if not constr_names == '':
                constrs_nfe_sample = constrs_run[:closest_nfe_idx,:]
                #constr_nfe_unique = constrs_nfe_sample[nfe_unique_idx,:]
                        
                if len(constr_runs_nfe) == 0:
                    #constr_runs_nfe = constr_nfe_unique
                    constr_runs_nfe = constrs_nfe_sample
                else:
                    #constr_runs_nfe = np.append(constr_runs_nfe, constr_nfe_unique, axis=0)
                    constr_runs_nfe = np.append(constr_runs_nfe, constrs_nfe_sample, axis=0)
        
        objs_runs_nfe_unique, current_unique_idx = np.unique(objs_runs_nfe, return_index=True, axis=0)
        if not constr_names == '':
            new_unique_idx = list(current_unique_idx.copy())
            for obj_idx in current_unique_idx:
                eq_obj_idx = [i for i in range(objs_runs_nfe.shape[0]) if (np.array_equal(objs_runs_nfe[i,:], objs_runs_nfe[obj_idx,:]) and i!=obj_idx)]
                if len(eq_obj_idx) > 0:
                    for current_eq_obj_idx in eq_obj_idx:
                        if all(constr_runs_nfe[current_eq_obj_idx,:] == 0) and any(constr_runs_nfe[obj_idx,:] != 0):
                            if obj_idx in new_unique_idx:
                                new_unique_idx.remove(obj_idx)
                            new_unique_idx.append(current_eq_obj_idx)
                            break # just need to add one instance of feasible design with same objectives
                      
            objs_runs_nfe_unique = objs_runs_nfe[new_unique_idx,:]
            constr_runs_nfe_unique = constr_runs_nfe[new_unique_idx,:]
            
            aggr_constrs_runs_nfe_unique = [np.mean(x) for x in constr_runs_nfe_unique]
            current_mooRunStats = MOORunStatistics(objs_runs_nfe_unique, aggr_constrs_runs_nfe_unique)
            pfs_idx_nfe = current_mooRunStats.compute_PF_idx_constrained(only_fullsat=True)
                
        else:
            
            current_mooRunStats = MOORunStatistics(objs_runs_nfe_unique)
            pfs_idx_nfe = current_mooRunStats.compute_PF_idx_unconstrained()
    
        nfe_pf_objs = objs_runs_nfe_unique[pfs_idx_nfe,:]
        
        nfe_sample_pfs['nfe: ' + str(nfe_sample)] = nfe_pf_objs
    
    case_pfs[case_key] = nfe_sample_pfs

# Plot Pareto Fronts for each case for each sample NFE (Assuming two objectives)
#obj_mult = [-1, 1] 
print('Plotting')
n_figs = len(nfe_samples)
for i in range(n_figs):
    fig = plt.figure()
    for case_key in list(case_objs.keys()):
        nfe_pfs = case_pfs[case_key]
        case_pf_nfe = nfe_pfs['nfe: ' + str(nfe_samples[i])]
                
        #pfs_obj1_case = [x[0] for x in case_pf_nfe]
        #pfs_obj2_case = [x[1] for x in case_pf_nfe]
        
        obj_mult = [1 for j in range(len(objs_max))]
        for k in range(len(obj_mult)): # if objective is negative (since its to be maximized), multiply by -1 for plotting
            if case_pf_nfe[0][k] < 0:
                obj_mult[k] = -1
            
        pfs_obj1 = [np.multiply(x[0], obj_mult[0]) for x in case_pf_nfe]
        pfs_obj2 = [np.multiply(x[1], obj_mult[1]) for x in case_pf_nfe]
        plt.scatter(pfs_obj1, pfs_obj2, c=plot_colors[case_key], marker=plot_markers[case_key], label=case_key)
    
    plt.xlabel(true_obj_names[0])
    plt.ylabel(true_obj_names[1])
    plt.legend(loc='best')
    plt.title('NFE = ' + str(nfe_samples[i]))
    
    