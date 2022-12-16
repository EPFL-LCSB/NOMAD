"""
Script to generate the plots for Figure 3 in the paper
- Plotting the mean growth, anthranilate, glucose for all the 10 models using the NRA suggested values when we allow:
- 2 fold changes in concentrations
- 3 fold changes in concentrations
- 10 fold changes in concentrations
- 20 fold changes in concentrations
- Unlimited changes (top yield control coefficients)
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np


plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})

CONCENTRATION_SCALING = 1e9
MAX_FOLD_CONC_CHANGES = [2, 3, 10, 20]

kinetic_models = ['1712,6', '1715,6', '1717,6', '2392,6', '2482,8', '3927,1', '4230,7', '4241,5', '4456,7',
                   '4468,6']
base_folder_to_data = '../../output/data/{}-fold-conc-change/'
folder_for_output = '../../output/figures/figure-3/'
if not os.path.exists(folder_for_output):
    os.makedirs(folder_for_output)

# Get the mean response values for the NRA designs for each allowable set of fold change in concentrations
list_mean_nra_designs = []
for this_fold_change in MAX_FOLD_CONC_CHANGES:
    df_nra_designs = []

    for this_model in kinetic_models:
        folder_to_data = base_folder_to_data.format(this_fold_change)
        filename = folder_to_data + '{}/ode_solutions.csv'.format(this_model)

        # Load dataframe of all solutions for the current kinetic model
        df = pd.read_csv(filename, index_col = 0)
        df_nra_designs.append(df[df['solution_id'] > 3])  # Only NRA design alternatives

    # Concatenate all the ODE solutions corresponding to the current set of allowable fold changes across all models
    df_nra_designs = pd.concat(df_nra_designs)
    list_mean_nra_designs.append(df_nra_designs.groupby('time').mean())

# Get the wildtype mean response and the MCA mean response from the 3-fold-conc-change folder
dfs_wt = []
dfs_mca = []
for this_model in kinetic_models:
    filename = '../../output/data/3-fold-conc-change/{}/ode_solutions.csv'.format(this_model)
    df = pd.read_csv(filename, index_col = 0)
    dfs_wt.append(df[df['solution_id'] == 0])  # Only wt
    dfs_mca.append(df[df['solution_id'] == 1])  # Only MCA

dfs_wt = pd.concat(dfs_wt)
dfs_mca = pd.concat(dfs_mca)
mean_wt = dfs_wt.groupby('time').mean()
mean_mca = dfs_mca.groupby('time').mean()

concentrations_to_plot = {'glc_D_e': 180.16 / CONCENTRATION_SCALING,
                          'anth_e': 136.13 / CONCENTRATION_SCALING,
                          'biomass_strain_1': 0.28e-12 /0.05 # Remember the reactor is 50 mL
                          }
ylabels = {'glc_D_e': 'Glucose (g/L)',
           'anth_e': 'Anthranilate (g/L)',
           'biomass_strain_1': 'Biomass (g/L)'
           }

markers = ['None', 'None', 'None', 'None', 'o']
linestyles = ['-', '--', '-.', ':', 'None']

for conc, scaling in concentrations_to_plot.items():

    # Plot wildtype for the current species concentration
    plt.plot(mean_wt.index, mean_wt[conc] * scaling,
             label='wt',
             color='orange', linewidth=2)

    # Plot the mean for each concentration constraint
    for ix, conc_change in enumerate(MAX_FOLD_CONC_CHANGES):
        df_to_plot = list_mean_nra_designs[ix]
        plt.plot(df_to_plot.index, df_to_plot[conc] * scaling,
                 label='{} fold change'.format(conc_change),
                 color='blue', linestyle=linestyles[ix], linewidth=2,
                 marker=markers[ix], markevery=50 + ix,  markersize=6, )

    # Plot mca for the current species concentration
    plt.plot(mean_mca.index, mean_mca[conc] * scaling,
             label='MCA',
             color='red', linewidth=2)

    plt.legend()
    plt.xlabel('Time (h)')
    plt.ylabel(ylabels[conc])
    plt.savefig(folder_for_output + 'Figure_3_{}.png'.format(conc))
    plt.close()

