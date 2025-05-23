"""
This is a script to generate the clustermap for the designs generated using the double mutant, eK_trpD9923_d2

 IMPORTANT!!
 - you need to install seaborn and sklearn to use this
 - type pip install scikit-learn
 - pip install seaborn
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn
from sklearn import *
import os


folder_to_data = './../../output/data/eK_trpD9923_d2/'
folder_for_output = './../../output/figures/figure-S.11/'
if not os.path.exists(folder_for_output):
    os.makedirs(folder_for_output)

# Read all the necesary files
df_designs = pd.read_csv(folder_to_data + 'all_unique_designs.csv', index_col = 0)
df_nra_sols = pd.read_csv(folder_to_data + 'all_nra_sols_across_models.csv', index_col=0)
df_subsystems = pd.read_csv('./../../data/general_rxn_subsystem.csv', index_col = 0)

# Get top designs by median / mean/ min - you can change this depending on what metric you want to use
indices_top_median = df_nra_sols.T.median().sort_values(ascending=False)[:5]
indices_top_mean = df_nra_sols.T.mean().sort_values(ascending=False)[:5]
indices_top_low = df_nra_sols.T.min().sort_values(ascending=False)[:5]

# Print top designs by mean value of solutions
for ix, value in indices_top_mean.iteritems():
    this_design = df_designs.iloc[ix]
    this_design = this_design[this_design == 1]
    print("{}, mean = {}".format(list(this_design.index), value))


'''
Clustermap plotting
'''
from matplotlib.colors import LinearSegmentedColormap

df_designs = df_designs.reset_index(drop=True)

# Get all the unique designs and their enzymes
df_designs = df_designs.dropna(axis=1, how='all')
df_designs = df_designs.fillna(0)

# Data reorganization - split into upregulations and downregulations
df_up = df_designs[[c for c in df_designs.columns if c.startswith('EU_')]]
df_down = df_designs[[c for c in df_designs.columns if c.startswith('ED_')]]
df_down = df_down * (-1)

# Combine up and down regulations into one df
df_up.columns = [c[3:] for c in df_up.columns]
df_down.columns = [c[3:] for c in df_down.columns]
col_list = list(set().union(df_up.columns, df_down.columns))
df_up = df_up.reindex(columns=col_list, fill_value=0)
df_down = df_down.reindex(columns=col_list, fill_value=0)

df_all = df_up + df_down
df_all = df_all.fillna(0)

# Assign and sort rows by subsystems
subs_assignment = [df_subsystems.loc[c][0] for c in df_all.columns]
df_all = df_all.reset_index(drop=True)
df_all_w_subs = df_all.T
df_all_w_subs['subsystem'] = subs_assignment
df_all_w_subs = df_all_w_subs.sort_values(by='subsystem', ascending = False)

# Colourmap for the rows based on subsystem
species = df_all_w_subs['subsystem']
colours = [ "red", "blue", "green", "yellow", "purple", "orange", "cyan","black","gray","magenta", "brown" ]
unique_subsystems = species.unique()
lut = dict(zip(unique_subsystems, colours))
row_colors = species.map(lut)

sns.set(font_scale = 3.5)
myColors = ((0.8, 0.0, 0.0, 1.0), (1.0, 1.0, 1.0, 1.0), (0.0, 0.0, 0.8, 1.0))
cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))
ax = sns.clustermap(df_all_w_subs.iloc[:,:-1], cmap=cmap, figsize = (30,35), row_cluster=False, row_colors = row_colors)
plt.tight_layout()
plt.savefig(folder_for_output + 'figure_S.11_enzyme_design.png')

