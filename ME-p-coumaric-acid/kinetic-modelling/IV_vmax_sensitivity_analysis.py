'''
Script to get stats on the ratio of Vmaxs for the Vmax sensitivity analysis.
We will analyse the distributions of the Vmaxs from the 2,000 models that are generated when we integrated
3x the exofluxomics from the fed-batch data. This uses information from the kinetic_parameter_sampling_vmax_sensitivty.py script.
But you can run this directly.
'''
from pytfa.io.json import load_json_model

from skimpy.io.yaml import load_yaml_model
from skimpy.analysis.oracle.load_pytfa_solution import load_fluxes, \
    load_concentrations, load_equilibrium_constants
from skimpy.core.parameters import ParameterValuePopulation
from skimpy.core.parameters import load_parameter_population
from skimpy.sampling.simple_parameter_sampler import SimpleParameterSampler
from skimpy.sampling.ga_parameter_sampler import GaParameterSampler
from skimpy.utils.general import get_stoichiometry
from skimpy.io.regulation import load_enzyme_regulation
from skimpy.analysis.mca.utils import get_dep_indep_vars_from_basis
import matplotlib.pyplot as plt

from sympy import Matrix
from scipy.sparse import csc_matrix as sparse_matrix

import pandas as pd
import numpy as np
from sys import argv
import glob


# Paths
path_to_old_param_set = './parameters/parameters_top_9_models.hdf5'
path_to_new_param_set = './parameters/vmax-sensitivity/parameters_sample_*.hdf5'

parameter_population = load_parameter_population(path_to_old_param_set)
old_param_set = parameter_population['3191,8'] # 3191,8 to keep the KMs
tmodel = load_json_model('./../models/constraint_based_model_scaffold.json')

# Get the Vmax info and ratios from all the generated parameter sets
dict_ = {}
for f in glob.iglob(path_to_new_param_set):
    ix_tfa = int(f.split('_sample_')[1][:-5])
    new_param_population = load_parameter_population(f)
    
    for ix_kin in new_param_population._index:
        # I want the ratio of Vmaxs compared to the orig parameter set
        dict_['{},{}'.format(ix_tfa, ix_kin)] = {k:v / old_param_set[k] for k, v in new_param_population[ix_kin].items() 
        if k.startswith('vmax_forward_')}
    
df_ratios = pd.DataFrame(dict_)
df_ratios_ = df_ratios.T

'''Print stats'''
stats = df_ratios_.agg('describe')
stats.to_csv('./output/vmax_ratio_stats_.csv')

stats_all = (df_ratios_).stack().agg('describe')
stats = stats.T

'''
Some stats on the number of reactions outside 75% and 125% of 3x the vmax for the perturbed models
'''
sample = pd.read_csv('./../data/steady_states.csv').iloc[3191] # Original steady-state profile
# Upper and lower bound 25% above and below
lb = 0.75
ub = 1.25

# 86 reactions had their median Vmaxs (across the 2,000 models) outside of the 25% range
outliers = stats[stats['50%'].le(lb) | stats['50%'].ge(ub)]
outlier_rxns = [i.replace('vmax_forward_', '') for i in outliers.index]

# 56 of these outliers were transport reactions
trans_rxns = [r for r in outlier_rxns if tmodel.reactions.get_by_id(r).thermo['isTrans']]

# The remaining 30 had a median Gamma (v_back / v_fwd) of 0.96 (very close to 1) meaning that they were close to equilibrium
nontrans_rxns = [r for r in outlier_rxns if r not in trans_rxns]
Gammas = np.exp(sample.loc[['LnGamma_' + i for i in nontrans_rxns]])
Gammas.agg('describe')
