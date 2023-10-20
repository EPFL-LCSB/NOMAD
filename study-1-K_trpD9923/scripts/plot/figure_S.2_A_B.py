"""
Plots for the nonlinear verification of the top 5 designs in a reactor setup
We plot the median/mean for wt and all the designs
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np


plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})

PERCENT_PERTURBATION = 25
CONCENTRATION_SCALING = 1e9
LIST_DESIGN_INDICES = [0, 1, 2, 3, 4]

# Choose folders according to which enzyme perturbation we want to plot for
folder_with_data = './../../output/verification-top-5-designs-reactor'
folder_for_output = './../../output/figures/figure-S.2/'
if not os.path.exists(folder_for_output):
    os.makedirs(folder_for_output)

# lists to store dataframes which will later be aggregated
dfs_wt_all = []
dfs_alt_all = []
time = np.linspace(0, 60, 1000)

# Get the mean response values for each design
list_mean_designs = []
for ix in range(0,5):

    filename = folder_with_data + '/design_{}_ode_sols.csv'.format(LIST_DESIGN_INDICES[ix])

    # Load dataframe of all solutions for the current design
    df = pd.read_csv(filename, index_col = 0)

    # Get the mean of the current design
    df_mean = df.groupby('time').mean()

    # Concat all the solutions for the top 5 designs to get the overall median
    dfs_alt_all.append(df)

    # Append the mean to the list
    list_mean_designs.append(df_mean)

# Get the median wildtype response
filename = './../../output/verification-top-5-designs-reactor/ode_sols_wt.csv'
df_wt = pd.read_csv(filename, index_col = 0)
mean_wt = df_wt.groupby('time').mean()

concentrations_to_plot = {'glc_D_e': 180.16 / CONCENTRATION_SCALING,
                          'anth_e': 136.13 / CONCENTRATION_SCALING,
                          'biomass_strain_1': 0.28e-12 /0.05,
                          }
ylabels = {'glc_D_e': 'Glucose (g/L)',
           'anth_e': 'Anthranilate (g/L)',
           'biomass_strain_1': 'Biomass (g/L)',
           }

markers = ['o', " ", " ", " ", " ", ]
line_styles = [ " ", '--', '-.', ':', '-']
for conc, scaling in concentrations_to_plot.items():

    # Plot wildtype for the current species
    plt.plot(mean_wt.index, mean_wt[conc] * scaling, label='wt', color='orange')

    # Plot the mean for each design
    for ix in range(5):
        df_to_plot = list_mean_designs[ix]
        plt.plot(df_to_plot.index, df_to_plot[conc] * scaling,
                 label = 'd-{}'.format(ix+1), color='blue',
                 marker=markers[ix], markevery=50, linestyle=line_styles[ix], linewidth=2)

    plt.legend()
    plt.xlabel('Time (h)')
    plt.ylabel(ylabels[conc])
    plt.savefig(folder_for_output + 'figure_S.2_A_B_{}.png'.format(conc))
    plt.close()

