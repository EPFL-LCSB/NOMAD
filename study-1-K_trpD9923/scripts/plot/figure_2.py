'''
Script to plot Figure 2:
- the responses of the 10 chosen kinetic models, K_trpD9923, that are representative of E.coli W3110 trpD9923
'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os


CONCENTRATION_SCALING = 1e9
kinetic_models = ['1712,6', '1715,6', '1717,6', '2392,6', '2482,8', '3927,1', '4230,7', '4241,5', '4456,7',
                   '4468,6']
folder_for_output = './../../output/figures/figure-2/'
if not os.path.exists(folder_for_output):
    os.makedirs(folder_for_output)

# Plotting parameters and data
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})
dict_scaling = {'biomass_strain_1' : 0.28e-12 / 0.05,
              'anth_e': 136.13 * 1e-9,
              'glc_D_e': 1e-9 * 180.156}
labels = {'biomass_strain_1' : 'Biomass (g)',
          'anth_e': 'Anthranilate (g/L)',
          'glc_D_e': 'Glucose (g/L)',
          }

# Load exp data
growth = pd.read_csv('./../../data/experiments_biomass.csv')
glc = pd.read_csv('./../../data/experiments_glc.csv')
anth = pd.read_csv('./../../data/experiments_anth.csv')
exp_data = {'biomass_strain_1' : growth,
            'anth_e': anth,
            'glc_D_e': glc,
            }

dfs_wt = []
# Collate all simulation results for the 10 models
for this_model in kinetic_models:
    filename = '../../output/data/3-fold-conc-change/{}/ode_solutions.csv'.format(this_model)
    df = pd.read_csv(filename, index_col = 0)
    dfs_wt.append(df[df['solution_id'] == 0])  # Only wt
dfs_wt = pd.concat(dfs_wt)

for conc, scaling in dict_scaling.items():

    # Get the exp data for the current concentration to plot
    exp_ = exp_data[conc]
    mean = exp_[exp_.columns[1]]
    time_exp = exp_[exp_.columns[0]]

    # Anth and glucose had error bounds
    if len(exp_.columns) > 2:
        lo = exp_[exp_.columns[2]]
        hi = exp_[exp_.columns[3]]
        plt.errorbar(time_exp, mean, yerr= np.asarray([mean-lo ,hi-mean]),
                     fmt='ko', label = 'exp', capsize= 5)
    else:
        plt.plot(time_exp, mean, 'ko', label = 'exp')

    # Plot simulated wildtype with the whole spread across the 10 models
    median_wt = dfs_wt.groupby('time').quantile(0.5)
    upper_wt = dfs_wt.groupby('time').quantile(0.99)
    lower_wt = dfs_wt.groupby('time').quantile(0.01)
    plt.plot(median_wt.index, median_wt[conc] * scaling, color='orange', label = 'sim')
    plt.fill_between(median_wt.index, lower_wt[conc] * scaling, upper_wt[conc] * scaling,
                     facecolor='orange',
                     interpolate=True,
                     alpha = 0.3)

    plt.legend()

    plt.xlabel('Time (h)')
    plt.ylabel(labels[conc])
    plt.xlim([0, 65])
    plt.savefig(folder_for_output + 'figure_2_{}.png'.format(conc))
    plt.close()

