"""
Plots for the nonlinear verification of the top 5 designs from eK_trpD9923_d2 in a reactor setup
We plot the median/mean for wt and all the designs
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import glob


plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})

PERCENT_PERTURBATION = 50
CONCENTRATION_SCALING = 1e9
LIST_DESIGN_INDICES = [0, 1, 2, 3, 4]

# Choose folders according to which enzyme perturbation we want to plot for
folder_with_data = './../../output/verification-top-5-designs-reactor/eK_trpD9923_d2'
folders_with_exp_strains = './../../output/data/eK_trpD9923_d2/**/'
path_to_design_sols = folder_with_data + '/design_*.csv'
folder_for_output = './../../output/figures/figure-S.13/'
if not os.path.exists(folder_for_output):
    os.makedirs(folder_for_output)

# Get mean curves for experimental strains using the different models
sols_wt = []
sols_d1 = []
sols_d2 = []
for this_folder in glob.iglob(folders_with_exp_strains):
    ode_sols = pd.read_csv(this_folder + 'ode_solutions.csv', header=0, index_col=0)
    sols_wt.append(ode_sols[ode_sols['solution_id'] == 0])
    sols_d1.append(ode_sols[ode_sols['solution_id'] == 2])
    sols_d2.append(ode_sols[ode_sols['solution_id'] == 3])

sols_wt = pd.concat(sols_wt, axis=0)
sols_d1 = pd.concat(sols_d1, axis=0)
sols_d2 = pd.concat(sols_d2, axis=0)

mean_wt = sols_wt.groupby('time').mean()
mean_d1 = sols_d1.groupby('time').mean()
mean_d2 = sols_d2.groupby('time').mean()

# lists to store dataframes which will later be aggregated
dfs_wt_all = []
dfs_alt_all = []
dfs_median = []
time = np.linspace(0, 60, 1000)

# Get the mean response values for each design
list_mean_designs = []
for filename in glob.iglob(path_to_design_sols):

    # filename = folder_with_data + '/design_{}_ode_sols.csv'.format(LIST_DESIGN_INDICES[ix])

    # Load dataframe of all solutions for the current design
    df = pd.read_csv(filename, index_col = 0)

    # Get the median of the current design
    df_mean = df.groupby('time').mean()

    # Append the median to the list
    list_mean_designs.append(df_mean)

# Get the median wildtype response
filename = folder_with_data + '/ode_sols_wt.csv'
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
source_data = []
for conc, scaling in concentrations_to_plot.items():

    # Plot wildtype for the current species
    plt.plot(mean_wt.index, mean_wt[conc] * scaling, label='wt', color='orange')

    # Plot d1 exp
    plt.plot(mean_d1.index, mean_d1[conc] * scaling, label='d1', color='red')

    # Plot d2 exp
    plt.plot(mean_d2.index, mean_d2[conc] * scaling, label='d2', color='black')
    df_temp = [mean_wt[conc] * scaling, mean_d1[conc] * scaling, mean_d2[conc] * scaling]

    # Plot the mean for each design
    for ix in range(5):
        df_to_plot = list_mean_designs[ix]
        plt.plot(df_to_plot.index, df_to_plot[conc] * scaling,
                 label = 'd-{}'.format(ix+1), color='blue',
                 marker=markers[ix], markevery=50, linestyle=line_styles[ix], linewidth=2)
        df_temp.append(df_to_plot[conc] * scaling)

    df_temp = pd.concat(df_temp, axis=1)
    df_temp.columns = [conc + i for i in ['_wt', '_ddpa', '_ddpa_tkt', '_d-1', '_d-2', '_d-3', '_d-4', '_d_5']]
    source_data.append(df_temp)

    # plt.legend()
    plt.xlabel('Time (h)')
    plt.ylabel(ylabels[conc])
    plt.savefig(folder_for_output + 'figure_S.13_a_{}.png'.format(conc))
    plt.close()

source_data = pd.concat(source_data, axis=1)
source_data.to_csv(folder_for_output + 'source_data_a.csv')

