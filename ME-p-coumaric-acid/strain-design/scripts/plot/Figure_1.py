"""
This is a script to generate the design vs enzyme clustermap and the NRA solution heatmap (i.e the NRA predicted
value for each combination of design and model) for Figure 1 in the manuscript

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


# Path to designs
path_to_designs = './../../output/designs/'

kinetic_models = ['3191,8', '3191,42', '3191,60', '3191,61', '3191,95', '3191,99',
                   '3191,104', '3191,139', '3191,142']

# Load all NRA designs for all models
list_designs = []
for this_model in kinetic_models:
    folder = path_to_designs + '{}'.format(this_model)
    df = pd.read_csv(folder + '/designs.csv', index_col = 0)
    list_designs.append(df)
df_designs = pd.concat(list_designs)

# Drop solutions and convert all designs to binary (1 if the enzyme is targeted in a design and 0 otherwise)
df_designs = df_designs.notnull().astype("int")
df_designs = df_designs.drop('solution', axis=1)

# Get only unique designs and convert into a list of designs
df_designs = df_designs.drop_duplicates()

'''
Create a clustermap with designs on the x axis and the enzymes that each design targets on the y axis
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

# Combine up and down regulations into one df - note that some enzymes have an up/down regulation depending on the design, these will be merged into a single enzyme
# So we might lose columns
df_up.columns = [c[3:] for c in df_up.columns]
df_down.columns = [c[3:] for c in df_down.columns]
col_list = list(set().union(df_up.columns, df_down.columns))
df_up = df_up.reindex(columns=col_list, fill_value=0)
df_down = df_down.reindex(columns=col_list, fill_value=0)

df_all = df_up + df_down
df_all = df_all.fillna(0)

# Order the columns for plotting
col_order = ['ME1m', 'CPR', 'TAL', 'ASPTA', 'DHQTi', 'CHORM', 'PPND2', 'DHQS', 'PGCD', 'PRFGS', 'GLNS', 'GLUDy',
             'GND', 'PGL', 'PRPPS', 'G6PDH2r', 'SUCOASm', 'ICDHym', 'ICDHyr', 'CSm', 'ICDHxm', 'MDH', 'MDHm', 'SUCD2_u6m',
             'PGK', 'PYRDC', 'AKGDam', 'PGM', 'PDHm', 'HEX1', 'PFK', 'r_4235', 'AKGDbm', 'PYK', 'PGI', 'FBA3',
             'ALDD2y', 'GCC2cm', 'PPA', 'NADH2_u6m']
df_all = df_all[col_order]

# Assign and sort rows by subsystems
df_all = df_all.reset_index(drop=True)

sns.set(font_scale = 3.5)
myColors = ((0.8, 0.0, 0.0, 1.0), (1.0, 1.0, 1.0, 1.0), (0.0, 0.0, 0.8, 1.0))
cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))
ax = sns.clustermap(df_all.T, cmap=cmap, figsize = (35,40), row_cluster=False)
plt.tight_layout()
plt.savefig('./../../output/figures/Figure_1.png', bbox_inches='tight')
plt.close()
