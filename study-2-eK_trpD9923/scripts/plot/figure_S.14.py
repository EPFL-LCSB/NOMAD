'''
This is a script to generate a figure with the old top 4 designs and the top 5 new 
NOMAD designs in one plot to show how close they are.
'''
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import glob


plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})
CONCENTRATION_SCALING = 1e9

path_to_new_data = './../../output/data/eK_trpD9923/**/'
path_to_new_top_designs = './../../output/verification-top-5-designs-reactor/eK_trpD9923/design_*.csv'
path_to_old_data = './../../output/data/old-designs-in-eK_trpD9923/design_{}_ode_sols.csv'
folder_for_output = '../../output/figures/figure-S.14/'
if not os.path.exists(folder_for_output):
    os.makedirs(folder_for_output)


'''
**** Load data from ALL new designs and get stats *****
'''
# Store data from new designs
dfs_wt_all, dfs_alt_all, dfs_ddpa_all, dfs_ddpa_tkt_all = [], [], [], []
for data_folder in glob.iglob(path_to_new_data):
    if not os.path.exists(data_folder):
        continue

    # Load dataframe of all solutions
    df = pd.read_csv(data_folder + 'ode_solutions.csv')
    df.columns = [c.replace('strain_1_','') for c in df.columns]
    n_sols = len(np.unique(df['solution_id']))

    # Append solutions according to type
    dfs_alt_all.append(df[df['solution_id']>3])     # Strain design alternatives
    dfs_wt_all.append(df[df['solution_id']==0])     # wild-type
    dfs_ddpa_all.append(df[df['solution_id'] == 2])  # Paper design - fbr DDPA
    dfs_ddpa_tkt_all.append(df[df['solution_id'] == 3])  # Paper design - fbr DDPA + TKT

dfs_wt_all = pd.concat(dfs_wt_all)
dfs_alt_all = pd.concat(dfs_alt_all)
dfs_ddpa_all = pd.concat(dfs_ddpa_all)
dfs_ddpa_tkt_all = pd.concat(dfs_ddpa_tkt_all)

def get_stats(df_):
    df_median = df_.groupby('time').median()
    df_upper = df_.groupby('time').quantile(0.75)
    df_lower = df_.groupby('time').quantile(0.25)
    return df_median, df_lower, df_upper

# Get mean, upper, and lower quartiles for each of the solution sets
dfs_wt_median, dfs_wt_lower, dfs_wt_upper = get_stats(dfs_wt_all)
dfs_alt_median, dfs_alt_lower, dfs_alt_upper = get_stats(dfs_alt_all)
dfs_ddpa_median, dfs_ddpa_lower, dfs_ddpa_upper = get_stats(dfs_ddpa_all)
dfs_ddpa_tkt_median, dfs_ddpa_tkt_lower, dfs_ddpa_tkt_upper = get_stats(dfs_ddpa_tkt_all)

'''
**** Load data from TOP FIVE new designs and get stats *****
'''
dfs_new_designs = []
for filename in glob.iglob(path_to_new_top_designs):
    df = pd.read_csv(filename)
    df.columns = [c.replace('strain_1_','') for c in df.columns]
    dfs_new_designs.append(df)

dfs_new_designs = pd.concat(dfs_new_designs)
dfs_new_median, dfs_new_lower, dfs_new_upper = get_stats(dfs_new_designs)

"""
************ Load data from top FOUR old designs and get stats *******************
"""
dfs_old_designs = []
for ix in range(4):
    df = pd.read_csv(path_to_old_data.format(ix))
    df.columns = [c.replace('strain_1_','') for c in df.columns]
    print("mean = {}".format(df.groupby('time').median().iloc[-1]['anth_e'] * 136.13/1e9))
    dfs_old_designs.append(df)

dfs_old_designs = pd.concat(dfs_old_designs)
dfs_old_median, dfs_old_lower, dfs_old_upper = get_stats(dfs_old_designs)

"""
**** PLOTTING ****
"""
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

    # Plot paper design - DDPA fbr + TKT
    plt.plot(time, dfs_ddpa_tkt_median[conc] * scaling, color='red', label='$aroG^{fbr}tktA$', linestyle='--',)

    # Plot old top 4 designs
    plt.plot(time, dfs_old_median[conc] * scaling, color='blue', label='top-4-old', linewidth=2)
    plt.fill_between(time, dfs_old_lower[conc] * scaling, dfs_old_upper[conc] * scaling, facecolor='blue',
                     interpolate=True,
                     alpha=0.1)

    # Plot new top 5 designs
    plt.plot(time, dfs_new_median[conc] * scaling, color='blue', label='top-5-new', linestyle='--')

    plt.legend(ncol=1)
    plt.xlabel('Time (h)')
    plt.xlim([0, 65])
    plt.ylabel(axis_labels[conc])
    plt.savefig(folder_for_output + 'figure_S.14_{}.png'.format(conc))
    plt.close()

