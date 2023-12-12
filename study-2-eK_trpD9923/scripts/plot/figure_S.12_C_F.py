"""
Plots for enzyme sensitivity analysis for top 5 designs from eK_trpD9923
- With all 3 enzymes being perturbed
- With each enzyme individually

Data is generated from V_expression_sensitivty_analysis.py. Report in Supplementary note IX
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import glob as glob


'''
Plotting options
- 0 = first enzyme (ANS), 
- 1 = second enzyme (DDPA) 
- 2 = third enzyme (ANS/PYK etc.)
- 3 for all 3 enzymes perturbed together
'''

plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})

PERCENT_PERTURBATION = 50
CONCENTRATION_SCALING = 1e9
ENZYME_NUMBERS_TO_PLOT = [2 ,3 ]

folder_for_output = './../../output/figures/figure-S.12/'
if not os.path.exists(folder_for_output):
    os.makedirs(folder_for_output)

figure_labels = ['d', 'e', 'f', 'c', ]
for enz_to_plot in ENZYME_NUMBERS_TO_PLOT:
    if enz_to_plot < 3:
        folder_with_data = './../../output/enzyme-sensitivity/eK_trpD9923/{}-percent/enzyme-{}/'.format(PERCENT_PERTURBATION, enz_to_plot)
    else:
        folder_with_data = './../../output/enzyme-sensitivity/eK_trpD9923/{}-percent/all-enzymes/'.format(PERCENT_PERTURBATION)
    path_to_data = folder_with_data + 'ode_sols_*.csv'
    # lists to store dataframes which will later be aggregated
    dfs_wt_all = []
    dfs_alt_all = []
    dfs_median = []
    time = np.linspace(0, 80, 1000)
    
    # Get the mean response values for each design - note that we leave out design 5 with ANS as the third enzymes
    # as it did not pass the test for robustness across models
    list_median_designs = []
    for filename in glob.iglob(path_to_data):

        # Load dataframe of all solutions for the current design
        df = pd.read_csv(filename, index_col = 0)

        # Get the mean of the current design
        df_mean = df.groupby('time').mean()

        # Append the median to the list
        list_median_designs.append(df_mean)

    # Get the median wildtype response
    filename = './../../output/verification-top-5-designs-reactor/eK_trpD9923/ode_sols_wt.csv'
    df_wt = pd.read_csv(filename, index_col = 0)

    median_wt = df_wt.groupby('time').mean()

    concentrations_to_plot = {'glc_D_e': 180.16 / CONCENTRATION_SCALING,
                              'anth_e': 136.13 / CONCENTRATION_SCALING,
                              'biomass_strain_1': 0.28e-12 /0.05
                              }
    ylabels = {'glc_D_e': 'Glucose (g/L)',
               'anth_e': 'Anthranilate (g/L)',
               'biomass_strain_1': 'Biomass (g/L)'
               }

    markers = ['o', " ", " ", " ", " ",]
    line_styles = [ " ", '--', '-.', ':', '-']
    source_data = []
    for conc, scaling in concentrations_to_plot.items():

        # Plot wildtype for the current species
        plt.plot(median_wt.index, median_wt[conc] * scaling, label='wt', color='orange')
        df_temp = [median_wt[conc] * scaling]

        # Plot the median for each design
        for ix in range(5):

            df_to_plot = list_median_designs[ix]
            plt.plot(df_to_plot.index, df_to_plot[conc] * scaling,
                     label = 'd-{}'.format(ix+1), color='blue',
                     marker=markers[ix], markevery=50, linestyle=line_styles[ix], linewidth=2)
            df_temp.append(df_to_plot[conc] * scaling)
            if conc == 'anth_e':
                print('Enz {}, design {} --> mean anthraniltae = {}'.format(enz_to_plot, ix, df_to_plot.iloc[-1][conc] * scaling))

        df_temp = pd.concat(df_temp, axis=1)
        df_temp.columns = [conc + i for i in ['_wt', '_d-1', '_d-2', '_d-3', '_d-4', '_d-5']]
        source_data.append(df_temp)

        # plt.legend()
        plt.xlabel('Time (h)')
        plt.xlim([0, 60])
        plt.ylabel(ylabels[conc])
        plt.savefig(folder_for_output + 'Figure_S.12_{}_{}.png'.format(figure_labels[enz_to_plot], conc))
        plt.close()

    source_data = pd.concat(source_data, axis=1)
    source_data.to_csv(folder_for_output + 'source_data_{}.csv'.format(figure_labels[enz_to_plot]))


