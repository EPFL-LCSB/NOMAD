'''
In this script, I plot the ODe solutions for the reactor simulations for the top 10 designs applied to the 9 kinetic models.
I only plot those which have PCA titer increases after 24 hours.

These are the plots used for Figures 2B, 2C and the Supplementary Figure S3
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


CONCENTRATION_SCALING = 1e9 # 1 mol to 1 mmol
DENSITY = 1200 # g/L
GDW_GWW_RATIO = 0.3 # Assumes 70% Water
CELL_VOLUME = 42 # mum3
REACTOR_VOLUME = 0.0025 #L
GDW_PER_CELL = GDW_GWW_RATIO * DENSITY * CELL_VOLUME * 1e-15 # 1e-15 L per mum3

FONT_SIZE = 20

# Paths - change them to reflect actual output folders
path_to_sols = './../../output/nonlinear-verification/'
folder_for_plots = './../../output/figures/'

os.makedirs(folder_for_plots, exist_ok=True)

concs_to_plot = ['glc__D_e', 'pca_e', 'biomass_strain_1']
labels = ['Glucose (g/L)', '$\it{p}$-coumaric acid (g/L)', 'Biomass (g/L)']
scaling = [180/CONCENTRATION_SCALING, 164/CONCENTRATION_SCALING, GDW_PER_CELL / REACTOR_VOLUME]

kinetic_models = ['3191,8', '3191,42', '3191,60', '3191,61', '3191,95', '3191,99',
                   '3191,104', '3191,139', '3191,142'] # 

# Load all the ODE simulation data
df_sols_wt = []
df_sols_designs = []
df_sols_designs_inc = []
for m in kinetic_models:
    try:
        sol_wt = pd.read_csv(path_to_sols + '{}/sol_wt.csv'.format(m), index_col=0)
        sol_designs = pd.read_csv(path_to_sols + '{}/sol_designs.csv'.format(m), index_col=0)
    except:
        print('Cannot find {}'.format(m))
        continue

    # some curation
    sol_designs['model'] = m
    sol_wt['time'] = sol_designs[sol_designs.solution_id == 0].time
    sol_wt['model'] = m
    
    pca_wt_24hr = sol_wt[(sol_wt.time - 24) == np.min(np.abs(sol_wt.time - 24))].pca_e.iloc[0]

    # Check if pca titer is higher at 24 hours
    pca_designs_24hr = sol_designs[(sol_designs.time - 24) == np.min(np.abs(sol_designs.time - 24))][['solution_id','pca_e']]
    ids_to_keep = pca_designs_24hr[pca_designs_24hr.pca_e.gt(pca_wt_24hr)].solution_id.values
    df_sols_wt.append(sol_wt)
    df_sols_designs.append(sol_designs)
    df_sols_designs_inc.append(sol_designs[sol_designs.solution_id.isin(ids_to_keep)])

df_sols_wt = pd.concat(df_sols_wt)
df_sols_designs = pd.concat(df_sols_designs)
df_sols_designs_inc = pd.concat(df_sols_designs_inc)

def get_quartiles(df):
    # First calculate median, lower and upper quartiles
    median = df.groupby('time').median()
    lower = df.groupby('time').quantile(0.25)
    upper = df.groupby('time').quantile(0.75)
    mean = df.groupby('time').mean()

    return mean, median, lower, upper

'''
Plot wildtype + all increasing strains 
'''

avg_wt, med_wt, lo_wt, up_wt  = get_quartiles(df_sols_wt)
avg_des_inc, med_des_inc, lo_des_inc, up_des_inc  = get_quartiles(df_sols_designs_inc)
for c, l, s in zip(concs_to_plot, labels, scaling):    

    # wt
    plt.plot(avg_wt.index, avg_wt[c] * s, color='orange')
    plt.fill_between(avg_wt.index, lo_wt[c] * s, up_wt[c] * s, facecolor='orange', label='ST10284',
                     interpolate=True,
                     alpha=0.1)

    # all designs
    plt.plot(avg_des_inc.index, avg_des_inc[c] * s, color='blue')
    plt.fill_between(avg_des_inc.index, lo_des_inc[c] * s, up_des_inc[c] * s, facecolor='blue', label='Designs',
                     interpolate=True,
                     alpha=0.1)

    plt.xlabel('Time (h)', fontsize=FONT_SIZE)
    plt.ylabel(l, fontsize=FONT_SIZE)
    plt.xticks(fontsize=12)  # Set all x-tick labels
    plt.yticks(fontsize=12)  # Set all x-tick labels

    plt.tight_layout()
    plt.savefig(folder_for_plots + 'Figure_2BC_{}.png'.format(c), bbox_inches='tight')
    plt.close()