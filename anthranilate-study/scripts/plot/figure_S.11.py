"""
Script to generate the plot for Supplementary Figure S.11 (Supplementary note VIII)
- Plot responses of the following
1. All wildtype solutions using these 10 models
2. Augmented in-silico version of aroG^fbr + tktA with 4 additional enzymes identified by NRA
3. Original in-silico version of aroG^fbr + tktA (5-fold upregulation)
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np


plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})

CONCENTRATION_SCALING = 1e9
MAX_FOLD_ENZ_ACTIVITY = 10 # OG was 5
MAX_FOLD_CONC_CHANGE = 5 # Changed from 3
MAX_N_ENZ = 6

kinetic_models = ['1712,6', '1715,6', '1717,6', '2392,6', '2482,8', '3927,1', '4230,7', '4241,5', '4456,7',
                   '4468,6']
folder_for_output = '../../output/figures/figure-S.11/'
if not os.path.exists(folder_for_output):
    os.makedirs(folder_for_output)

# lists to store dataframes which will later be aggregated
dfs_wt_all = []
dfs_ddpa = []
dfs_ddpa_tkt = []
list_dfs_ddpa_tkt_shiki = []

for this_model in kinetic_models:
    # data_folder = './../../output/rebuttal-data/exp-val-5fold/{}'.format(this_model)
    data_folder = './../../output/data/S.VIII-aroGtktA-analysis/reactor-verification/{}'.format(this_model)

    if not os.path.exists(data_folder):
        continue

    # Load dataframe of all solutions
    df = pd.read_csv(data_folder+'/ode_solutions.csv')
    df.columns = [c.replace('strain_1_','') for c in df.columns]
    n_sols = len(np.unique(df['solution_id']))

    # Append solutions according to type
    dfs_wt_all.append(df[df['solution_id']==0])     # wild-type
    dfs_ddpa.append(df[df['solution_id'] == 1])  # Paper design - fbr DDPA
    dfs_ddpa_tkt.append(df[df['solution_id'] == 2])  # Paper design - fbr DDPA
    list_dfs_ddpa_tkt_shiki.append(df[df['solution_id'] == 3])  # Paper design - fbr DDPA + TKT

# Get mean, upper, and lower quartiles for each of the solution sets
dfs_wt_all = pd.concat(dfs_wt_all)
dfs_wt_median = dfs_wt_all.groupby('time').quantile(0.5)
dfs_wt_mean = dfs_wt_all.groupby('time').mean()
dfs_wt_upper = dfs_wt_all.groupby('time').quantile(0.75)
dfs_wt_lower = dfs_wt_all.groupby('time').quantile(0.25)

dfs_ddpa = pd.concat(dfs_ddpa)
dfs_ddpa_median = dfs_ddpa.groupby('time').quantile(0.5)
dfs_ddpa_upper = dfs_ddpa.groupby('time').quantile(0.75)
dfs_ddpa_lower = dfs_ddpa.groupby('time').quantile(0.25)

dfs_ddpa_tkt = pd.concat(dfs_ddpa_tkt)
dfs_ddpa_tkt_median = dfs_ddpa_tkt.groupby('time').quantile(0.5)
dfs_ddpa_tkt_mean = dfs_ddpa_tkt.groupby('time').mean()
dfs_ddpa_tkt_upper = dfs_ddpa_tkt.groupby('time').quantile(0.75)
dfs_ddpa_tkt_lower = dfs_ddpa_tkt.groupby('time').quantile(0.25)

dfs_ddpa_tkt_shiki = pd.concat(list_dfs_ddpa_tkt_shiki)
dfs_ddpa_tkt_shiki_median = dfs_ddpa_tkt_shiki.groupby('time').quantile(0.5)
dfs_ddpa_tkt_shiki_mean = dfs_ddpa_tkt_shiki.groupby('time').mean()
dfs_ddpa_tkt_shiki_upper = dfs_ddpa_tkt_shiki.groupby('time').quantile(0.75)
dfs_ddpa_tkt_shiki_lower = dfs_ddpa_tkt_shiki.groupby('time').quantile(0.25)


concentrations_to_plot = {'glc_D_e': 180.16 / CONCENTRATION_SCALING,
                          'anth_e': 136.13 / CONCENTRATION_SCALING,
                          'biomass_strain_1': 0.28e-12 / 0.05,
                          }
# concentrations_to_plot = {'anth_e': 136.13 / CONCENTRATION_SCALING,
#                           }


axis_labels = {'glc_D_e': 'Glucose (g/L)',
               'anth_e': 'Anthranilate (g/L)',
               'biomass_strain_1': 'Biomass (g/L)',
               }


time = list(dfs_wt_median.index)
for conc, scaling in concentrations_to_plot.items():

    # Plot wildtype
    plt.plot(time, dfs_wt_median[conc] * scaling, color='orange', label = 'wt')
    # plt.fill_between(time, dfs_wt_lower[conc] * scaling, dfs_wt_upper[conc] * scaling, facecolor='orange', interpolate=True,
    #                  alpha = 0.3)

    # Plot design alternatives
    plt.plot(time, dfs_ddpa_tkt_shiki_median[conc] * scaling, color='blue', label = '$aroG^{fbr}tktA+$')
    plt.fill_between(time, dfs_ddpa_tkt_shiki_lower[conc] * scaling, dfs_ddpa_tkt_shiki_upper[conc] * scaling, facecolor='blue', interpolate=True,
                     alpha=0.1)

    # Plot MCA based approach
    # plt.plot(time, dfs_fcc_median[conc] * scaling, color='red', label='mca')
    # plt.fill_between(time, dfs_fcc_lower[conc] * scaling, dfs_fcc_upper[conc] * scaling, facecolor='red',
    #                  interpolate=True,
    #                  alpha=0.1)

    # Plot paper design - DDPA fbr
    # plt.plot(time, dfs_ddpa_median[conc] * scaling, color='black', label='$aroG^{fbr}$', linestyle='--',)
    # plt.fill_between(time, dfs_ddpa_lower[conc] * scaling, dfs_ddpa_upper[conc] * scaling, facecolor='black',
    #                  interpolate=True,
    #                  alpha=0.1)

    # Plot paper design - DDPA fbr + TKT
    plt.plot(time, dfs_ddpa_tkt_median[conc] * scaling, color='red', label='$aroG^{fbr}tktA-OLD$', linestyle='--',)
    # plt.fill_between(time, dfs_ddpa_tkt_lower[conc] * scaling, dfs_ddpa_tkt_upper[conc] * scaling, facecolor='red',
    #                  interpolate=True,
    #                  alpha=0.1)

    plt.legend()
    plt.xlabel('Time (h)')
    plt.xlim([0, 65])
    plt.ylabel(axis_labels[conc])
    plt.savefig(folder_for_output + '/Figure_S.11_{}.png'.format(conc))
    plt.close()

# Plot individual curves
#
# time = list(dfs_wt_median.index)
# for conc, scaling in concentrations_to_plot.items():
#
#
#     # Plot wildtype
#     plt.plot(time, dfs_wt_median[conc] * scaling, color='orange', label = 'wt')
#     # plt.fill_between(time, dfs_wt_lower[conc] * scaling, dfs_wt_upper[conc] * scaling, facecolor='orange', interpolate=True,
#     #                  alpha = 0.3)
#
#     for i_, this_sol in enumerate(list_dfs_ddpa_tkt_shiki):
#         # Plot design alternatives
#         if i_ == 9:
#             plt.plot(time, this_sol[conc] * scaling, color='blue', label = '$aroG^{fbr}tktA+$')
#         else:
#             plt.plot(time, this_sol[conc] * scaling, color='blue')
#
#     # Plot paper design - DDPA fbr + TKT
#     plt.plot(time, dfs_ddpa_tkt_median[conc] * scaling, color='red', label='$aroG^{fbr}tktA-OLD$', linestyle='--',)
#
#     plt.legend()
#     plt.xlabel('Time (h)')
#     plt.xlim([0, 65])
#     plt.ylabel(axis_labels[conc])
#     plt.savefig(folder_for_output + '/Figure_5_{}_shiki_w_ddpa_indi.png'.format(conc))
#     plt.close()