'''
- This is a script that collates all the designs that have been generated for ddpa + tkt + 4 more enzymes to address
the reviewer's comments
- Output is the top design along with the average fold changes for each of the enzymes in the top design
- NOTE: you don't need to run this script. I used it just to ascertain the mean fold changes for the top design.
- The top design targets DDPA, TKT1, TKT2, DHQS, SHKK, ANS, and GLUDy
'''
import pandas as pd
import numpy as np


# Paths
path_to_data = './../output/data/S.VIII-aroGtktA-analysis/design-generation/{}'
kinetic_models = ['1712,6', '1715,6', '1717,6', '2392,6', '2482,8', '3927,1', '4230,7', '4241,5', '4456,7', '4468,6']

# Collate all the designs into one dataframe with all the designs and the corresponding enz activity changes and NRA solution
designs = []
for m_ in kinetic_models:
    des_ = pd.read_csv(path_to_data.format(m_) + '/designs.csv', header=0, index_col=0)
    designs.append(des_)
df_designs = pd.concat(designs)

# Convert the dataframe to a binary one - 1 if the enzyme is in the design (non-zero enz activity change) and 0 otherwise
df_designs_binary = df_designs.copy()
df_designs_binary[~df_designs_binary.isnull()] = 1
df_designs_binary[df_designs_binary.isnull()] = 0

# Get count duplicates multiple columns using dataframe.pivot_table()
# This shows us which is the most frequently occuring design (the last one)
# Design with the top recurrence is EU_ANS, EU_DDPA, EU_DHQS, EU_SHKK, EU_TKT1, EU_TKT2, ED_GLUDy
df_recurrence = df_designs_binary.pivot_table(index = list(df_designs_binary.columns), aggfunc ='size')
print(df_recurrence)

# Get average values of fold changes for this top design across the 6 models for which it appears
top_design_cols = ['EU_ANS', 'EU_DDPA', 'EU_DHQS', 'EU_SHKK', 'EU_TKT1', 'EU_TKT2', 'ED_GLUDy']
df_top_design = df_designs[top_design_cols].dropna()
print(np.exp(df_top_design.mean()))




