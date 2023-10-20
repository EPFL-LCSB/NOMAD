'''
Script to plot Figure 6B:
- the responses of the 35 models calibrated on the behavior of trpD9923 and trpD9923/pJLaroGfbr
'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os


CONCENTRATION_SCALING = 1e9
styles = ['ko', 'kx', 'k^']
folder_to_exp_data = './../../data/fermentation-experiments/'
path_to_wt_sols = './../../output/data/top-35-models/ode_sols_wt.csv'
path_to_d1_sols = './../../output/data/top-35-models/ode_sols_d1.csv'
path_to_d2_sols = './../../output/data/top-35-models/ode_sols_d2.csv'
folder_for_output = './../../output/figures/figure-6/'
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
growth_wt = pd.read_csv(folder_to_exp_data + 'biomass_wt.csv')
growth_d1 = pd.read_csv(folder_to_exp_data + 'biomass_aroG.csv')
growth_d2 = pd.read_csv(folder_to_exp_data + 'biomass_aroG_tktA.csv')
glc_wt = pd.read_csv(folder_to_exp_data + 'glc_wt.csv')
glc_d1 = pd.read_csv(folder_to_exp_data + 'glc_aroG.csv')
glc_d2 = pd.read_csv(folder_to_exp_data + 'glc_aroG_tktA.csv')
anth_wt = pd.read_csv(folder_to_exp_data + 'anth_wt.csv')
anth_d1 = pd.read_csv(folder_to_exp_data + 'anth_aroG_curated.csv')
anth_d2 = pd.read_csv(folder_to_exp_data + 'anth_aroG_tktA_curated.csv')

exp_data = {'biomass_strain_1' : [growth_wt, growth_d1, growth_d2],
              'anth_e': [anth_wt, anth_d1, anth_d2],
              'glc_D_e': [glc_wt, glc_d1, glc_d2]}

# Get simulation results
dfs_wt = pd.read_csv(path_to_wt_sols, header=0, index_col=0)
dfs_d1 = pd.read_csv(path_to_d1_sols, header=0, index_col=0)
dfs_d2 = pd.read_csv(path_to_d2_sols, header=0, index_col=0)

# Calculate median and lower and upper quartile of the responses
median_wt = dfs_wt.groupby('time').quantile(0.5)
upper_wt = dfs_wt.groupby('time').quantile(0.75)
lower_wt = dfs_wt.groupby('time').quantile(0.25)
median_d1 = dfs_d1.groupby('time').quantile(0.5)
upper_d1 = dfs_d1.groupby('time').quantile(0.75)
lower_d1 = dfs_d1.groupby('time').quantile(0.25)
median_d2 = dfs_d2.groupby('time').quantile(0.5)
upper_d2 = dfs_d2.groupby('time').quantile(0.75)
lower_d2 = dfs_d2.groupby('time').quantile(0.25)

# ******************************************************************************************
for conc, scaling in dict_scaling.items():
    
    # Plot the experimental data first for each strain
    for i, data_ in enumerate(exp_data[conc]):
        mean = data_[data_.columns[1]]
        time_exp = data_[data_.columns[0]]

        # Some data have error bars
        if len(data_.columns) > 2:
            lo = data_[data_.columns[2]]
            hi = data_[data_.columns[3]]
            plt.errorbar(time_exp, mean, yerr= np.asarray([mean-lo ,hi-mean]),
                        fmt=styles[i], label = 'exp', capsize= 5)
        else:
            plt.plot(time_exp, mean, styles[i], label = 'exp')

    # Plot trpD9923^sim-enh
    plt.plot(median_wt.index, median_wt[conc] * scaling, color='orange', label = 'sim')
    plt.fill_between(median_wt.index, lower_wt[conc] * scaling, upper_wt[conc] * scaling,
                     facecolor='orange',
                     interpolate=True,
                     alpha = 0.2)

    # Plot trpD9923^sim-enh/pJLaroG^fbr 
    plt.plot(median_d1.index, median_d1[conc] * scaling, color='black', label = 'sim')
    plt.fill_between(median_d1.index, lower_d1[conc] * scaling, upper_d1[conc] * scaling,
                     facecolor='black',
                     interpolate=True,
                     alpha = 0.1)
    
    # Plot trpD9923^sim-enh/pJLaroG^fbrtktA
    plt.plot(median_d2.index, median_d2[conc] * scaling, color='red', label = 'sim')
    plt.fill_between(median_d2.index, lower_d2[conc] * scaling, upper_d2[conc] * scaling,
                     facecolor='red',
                     interpolate=True,
                     alpha = 0.1)
    
    # plt.legend()
    plt.xlabel('Time (h)')
    plt.xlim([0, 65])
    plt.ylabel(labels[conc])
    plt.savefig(folder_for_output + '/figure_6A_{}.png'.format(conc))
    plt.close()

