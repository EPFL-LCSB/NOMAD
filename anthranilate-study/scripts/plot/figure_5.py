"""
Script to generate the plot for Figure 5
- Plot responses of the following
1. All NOMAD generated solutions across all the 10 models
2. All wildtype solutions usnig these 10 models
3. Paper design 1 with a feedback resistant version of aroG
4. Paper design 2 with a feedback resistant version of aroG + upregulation of transketolase
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
folder_for_output = '../../output/figures/figure-5/'
if not os.path.exists(folder_for_output):
    os.makedirs(folder_for_output)

# lists to store dataframes which will later be aggregated
dfs_wt_all = []
dfs_alt_all = []
dfs_fcc_all = []
dfs_ddpa_all = []
dfs_ddpa_tkt_all = []

for this_model in kinetic_models:
    data_folder = './../../output/data/3-fold-conc-change/{}'.format(this_model)
    if not os.path.exists(data_folder):
        continue

    # Load dataframe of all solutions
    df = pd.read_csv(data_folder+'/ode_solutions.csv')
    df.columns = [c.replace('strain_1_','') for c in df.columns]
    n_sols = len(np.unique(df['solution_id']))

    # Append solutions according to type
    dfs_alt_all.append(df[df['solution_id']>3])     # Strain design alternatives
    dfs_wt_all.append(df[df['solution_id']==0])     # wild-type
    dfs_fcc_all.append(df[df['solution_id'] == 1])  # Using top yield CCs
    dfs_ddpa_all.append(df[df['solution_id'] == 2])  # Paper design - fbr DDPA
    dfs_ddpa_tkt_all.append(df[df['solution_id'] == 3])  # Paper design - fbr DDPA + TKT

# Get mean, upper, and lower quartiles for each of the solution sets
dfs_wt_all = pd.concat(dfs_wt_all)
dfs_wt_median = dfs_wt_all.groupby('time').quantile(0.5)
dfs_wt_upper = dfs_wt_all.groupby('time').quantile(0.75)
dfs_wt_lower = dfs_wt_all.groupby('time').quantile(0.25)

dfs_alt_all = pd.concat(dfs_alt_all)
dfs_alt_median = dfs_alt_all.groupby('time').quantile(0.5)
dfs_alt_upper = dfs_alt_all.groupby('time').quantile(0.75)
dfs_alt_lower = dfs_alt_all.groupby('time').quantile(0.25)

dfs_fcc_all = pd.concat(dfs_fcc_all)
dfs_fcc_median = dfs_fcc_all.groupby('time').quantile(0.5)
dfs_fcc_upper = dfs_fcc_all.groupby('time').quantile(0.75)
dfs_fcc_lower = dfs_fcc_all.groupby('time').quantile(0.25)

dfs_ddpa_all = pd.concat(dfs_ddpa_all)
dfs_ddpa_median = dfs_ddpa_all.groupby('time').quantile(0.5)
dfs_ddpa_upper = dfs_ddpa_all.groupby('time').quantile(0.75)
dfs_ddpa_lower = dfs_ddpa_all.groupby('time').quantile(0.25)

dfs_ddpa_tkt_all = pd.concat(dfs_ddpa_tkt_all)
dfs_ddpa_tkt_median = dfs_ddpa_tkt_all.groupby('time').quantile(0.5)
dfs_ddpa_tkt_upper = dfs_ddpa_tkt_all.groupby('time').quantile(0.75)
dfs_ddpa_tkt_lower = dfs_ddpa_tkt_all.groupby('time').quantile(0.25)

concentrations_to_plot = {'glc_D_e': 180.16 / CONCENTRATION_SCALING,
                          'anth_e': 136.13 / CONCENTRATION_SCALING,
                          'biomass_strain_1': 0.28e-12 / 0.05,
                          }

axis_labels = {'glc_D_e': 'Glucose (g/L)',
               'anth_e': 'Anthranilate (g/L)',
               'biomass_strain_1': 'Biomass (g/L)',
               }


time = list(dfs_wt_median.index)
for conc, scaling in concentrations_to_plot.items():

    # Plot wildtype
    plt.plot(time, dfs_wt_median[conc] * scaling, color='orange', label = 'wt')
    plt.fill_between(time, dfs_wt_lower[conc] * scaling, dfs_wt_upper[conc] * scaling, facecolor='orange', interpolate=True,
                     alpha = 0.3)

    # Plot design alternatives
    plt.plot(time, dfs_alt_median[conc] * scaling, color='blue', label = 'NOMAD')
    plt.fill_between(time, dfs_alt_lower[conc] * scaling, dfs_alt_upper[conc] * scaling, facecolor='blue', interpolate=True,
                     alpha=0.1)

    # Plot MCA based approach
    # plt.plot(time, dfs_fcc_median[conc] * scaling, color='red', label='mca')
    # plt.fill_between(time, dfs_fcc_lower[conc] * scaling, dfs_fcc_upper[conc] * scaling, facecolor='red',
    #                  interpolate=True,
    #                  alpha=0.1)

    # Plot paper design - DDPA fbr
    plt.plot(time, dfs_ddpa_median[conc] * scaling, color='black', label='$aroG^{fbr}$', linestyle='--',)
    plt.fill_between(time, dfs_ddpa_lower[conc] * scaling, dfs_ddpa_upper[conc] * scaling, facecolor='black',
                     interpolate=True,
                     alpha=0.1)

    # Plot paper design - DDPA fbr + TKT
    plt.plot(time, dfs_ddpa_tkt_median[conc] * scaling, color='red', label='$aroG^{fbr}tktA$', linestyle='--',)
    plt.fill_between(time, dfs_ddpa_tkt_lower[conc] * scaling, dfs_ddpa_tkt_upper[conc] * scaling, facecolor='red',
                     interpolate=True,
                     alpha=0.1)

    plt.legend()
    plt.xlabel('Time (h)')
    plt.xlim([0, 65])
    plt.ylabel(axis_labels[conc])
    plt.savefig(folder_for_output + 'Figure_5_{}.png'.format(conc))
    plt.close()

