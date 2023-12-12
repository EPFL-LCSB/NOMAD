"""
Script to generate the plots for Supplementary Note II, Figures S.2 to S.4
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np


plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})

CONCENTRATION_SCALING = 1e9
kinetic_models = ['1712,6', '1715,6', '1717,6', '2392,6', '2482,8', '3927,1', '4230,7', '4241,5', '4456,7',
                   '4468,6']
folder_for_output = '../../output/figures/figure-S.2_S.4/'
if not os.path.exists(folder_for_output):
    os.makedirs(folder_for_output)

# lists to store dataframes which will later be aggregated
dfs_wt_all = []
dfs_mca_1_2 = []
dfs_mca_2 = []
dfs_mca_5 = []

for this_model in kinetic_models:
    data_folder = './../../output/mca-study/{}'.format(this_model)
    if not os.path.exists(data_folder):
        continue

    # Load dataframe of all solutions
    df = pd.read_csv(data_folder + '/ode_solutions.csv')
    df.columns = [c.replace('strain_1_','') for c in df.columns]
    n_sols = len(np.unique(df['solution_id']))

    # Append solutions according to type
    dfs_wt_all.append(df[df['solution_id']==0])     # wild-type
    dfs_mca_1_2.append(df[df['solution_id'] == 1])  # Using top yield CCs w 1.2 fold changes in enzyme activity
    dfs_mca_2.append(df[df['solution_id'] == 2])  # Using top yield CCs w 2 fold changes
    dfs_mca_5.append(df[df['solution_id'] == 3])  # Using top yield CCs w 5 fold changes

# Get mean, upper, and lower quartiles for the wt
dfs_wt_all = pd.concat(dfs_wt_all)
dfs_wt_median = dfs_wt_all.groupby('time').quantile(0.5)
dfs_wt_upper = dfs_wt_all.groupby('time').quantile(0.75)
dfs_wt_lower = dfs_wt_all.groupby('time').quantile(0.25)
time = list(dfs_wt_median.index)

dfs_mca_1_2 = pd.concat(dfs_mca_1_2)
dfs_mca_2 = pd.concat(dfs_mca_2)
dfs_mca_5 = pd.concat(dfs_mca_5)

concentrations_to_plot = {'glc_D_e': 180.16 / CONCENTRATION_SCALING,
                          'anth_e': 136.13 / CONCENTRATION_SCALING,
                          'biomass_strain_1': 0.28e-12 / 0.05,
                          }

axis_labels = {'glc_D_e': 'Glucose (g/L)',
               'anth_e': 'Anthranilate (g/L)',
               'biomass_strain_1': 'Biomass (g/L)',
               }
ENZ_ACTIVITIES = ['1.2', '2', '5']
figure_names = ['S.2', 'S.3', 'S.4']
source_data = []
for this_enz_act, figure_name, this_df in zip(ENZ_ACTIVITIES, figure_names, [dfs_mca_1_2, dfs_mca_2, dfs_mca_5]):

    df_median = this_df.groupby('time').quantile(0.5)
    df_upper = this_df.groupby('time').quantile(0.75)
    df_lower = this_df.groupby('time').quantile(0.25)

    for conc, scaling in concentrations_to_plot.items():

        # Plot wildtype
        plt.plot(time, dfs_wt_median[conc] * scaling, color='orange', label = 'wt')
        plt.fill_between(time, dfs_wt_lower[conc] * scaling, dfs_wt_upper[conc] * scaling, facecolor='orange', interpolate=True,
                         alpha = 0.3)

        # Plot design alternatives
        plt.plot(time, df_median[conc] * scaling, color='red', label = 'mca {}-fold'.format(this_enz_act))
        plt.fill_between(time, df_lower[conc] * scaling, df_upper[conc] * scaling, facecolor='red', interpolate=True, alpha=0.1)

        df_temp = pd.concat([dfs_wt_median[conc] * scaling, dfs_wt_lower[conc] * scaling, dfs_wt_upper[conc] * scaling,
                             df_median[conc] * scaling, df_lower[conc] * scaling, df_upper[conc] * scaling,
                             ],
                            axis=1)
        df_temp.columns = [conc + '_' + i + '_{}_fold_change'.format(this_enz_act)
                           for i in ['wt_median', 'wt_lower', 'wt_upper', 'designs_median', 'designs_lower', 'designs_upper']
                           ]
        source_data.append(df_temp)
        plt.legend()
        plt.xlabel('Time (h)')
        plt.xlim([0, 65])
        plt.ylabel(axis_labels[conc])
        plt.savefig(folder_for_output + 'figure_{}_{}.png'.format(figure_name, conc))
        plt.close()

source_data = pd.concat(source_data, axis=1)
source_data.to_csv(folder_for_output + 'source_data.csv')