'''
Script to plot the figure S.8 - performance of the double mutant targeting DDPA and GLUDy only
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


# Paths
folder_to_data = './../../output/data/S.VII-double-mutant'
folder_for_output = './../../output/figures/figure-S.8'

if not os.path.exists(folder_for_output):
    os.makedirs(folder_for_output)

# Load data
df_sols_wt = pd.read_csv(folder_to_data + '/ode_sols_wt.csv', header=0, index_col=0)
df_sols_design = pd.read_csv(folder_to_data + '/ode_sols_ddpa_gludy.csv', header=0, index_col=0)

df_wt_mean = df_sols_wt.groupby('time').quantile(0.5)
df_wt_upper = df_sols_wt.groupby('time').quantile(0.75)
df_wt_lower = df_sols_wt.groupby('time').quantile(0.25)

df_designs_mean = df_sols_design.groupby('time').quantile(0.5)
df_designs_upper = df_sols_design.groupby('time').quantile(0.75)
df_designs_lower = df_sols_design.groupby('time').quantile(0.25)

concentrations_to_plot = {'glc_D_e': 1 / 1e9,
                          'anth_e': 136.13 / 1e9,
                          'biomass_strain_1': 0.28e-12 / 0.05}

ylabels = {'glc_D_e': 'Glucose (g/L)',
          'anth_e': 'Anthranilate (g/L)',
          'biomass_strain_1': 'Biomass (g/L)',
          }

for conc, scaling in concentrations_to_plot.items():
    time = df_wt_mean.index
    # Plot wildtype
    plt.plot(time, df_wt_mean[conc] * scaling, color='orange', label='wt')
    plt.fill_between(time, df_wt_lower[conc] * scaling, df_wt_upper[conc] * scaling, facecolor='orange',
                     interpolate=True,
                     alpha=0.1)

    # Plot design alternatives
    plt.plot(time, df_designs_mean[conc] * scaling, color='blue', label='NOMAD - DDPA + GLUDy')
    plt.fill_between(time, df_designs_lower[conc] * scaling, df_designs_upper[conc] * scaling,
                     facecolor='blue',
                     interpolate=True,
                     alpha=0.1)

    plt.legend()
    plt.xlabel('Time (h)')
    plt.ylabel(ylabels[conc])
    plt.savefig(folder_for_output + '/figure_S.8_{}.png'.format(conc))
    plt.close()
