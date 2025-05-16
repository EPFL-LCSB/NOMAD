'''
Script to plot a histogram of Vmaxs and the ratio for the Vmax sensitivity analysis (Supplementary Figures)
This is for Supplementary Figure S2
'''
from pytfa.io.json import load_json_model
from skimpy.core.parameters import load_parameter_population
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import glob


# Paths
path_to_old_param_set = './../../../kinetic-modelling/parameters/parameters_top_9_models.hdf5'
path_to_new_param_set = './../../../kinetic-modelling/parameters/vmax-sensitivity/parameters_sample_*.hdf5'
path_to_tmodel = './../../../models/constraint_based_model_scaffold.json'

parameter_population = load_parameter_population(path_to_old_param_set)
old_param_set = parameter_population['3191,8'] # 3191,8 to keep the KMs
tmodel = load_json_model(path_to_tmodel)

# Get the Vmax info and ratios from all the generated parameter sets
dict_ = {}
for f in glob.iglob(path_to_new_param_set):
    ix_tfa = int(f.split('_sample_')[1][:-5])
    new_param_population = load_parameter_population(f)
    
    for ix_kin in new_param_population._index:
        # I want the ratio of Vmaxs compared to the orig parameter set without 3x fluxes
        dict_['{},{}'.format(ix_tfa, ix_kin)] = {k:v / old_param_set[k] for k, v in new_param_population[ix_kin].items() 
        if k.startswith('vmax_forward_')}
    
df_ratios = pd.DataFrame(dict_)
df_ratios_ = df_ratios.T

''' Plot histogram '''
stats_all = (df_ratios_).stack().agg('describe')
stats = stats.T

# Upper and lower bound 25% above and below
lb = 0.75
ub = 1.25

plt.hist(np.log(df_ratios_).stack(), bins=300, alpha=0.8)  # alpha controls the transparency
plt.plot([np.log(lb), np.log(lb)], [0, 70000], 'k-')
plt.plot([np.log(ub), np.log(ub)], [0, 70000], 'k-')
plt.xlabel('$V_{max}^k / (3V_{max}^{fb})$', fontsize=16)
plt.ylabel('Frequency', fontsize=16)
plt.tight_layout()
plt.savefig('./../../output/figures/Figure_S2.png', bbox_inches='tight')
plt.close()
