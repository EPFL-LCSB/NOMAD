'''
This is a script that plots subplot A of Figure 3 in the main paper
'''
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


df = pd.read_csv('./../../output/nonlinear-verification/performance_analysis.csv', index_col=0)

# This is for plotting purposes - a link between the in-silico design indices and their experimental strain counterparts
comp_exp_links = {  0: 'ST14204',
                    1: 'ST14210',
                    4: 'ST14206',
                    5: 'ST14201',
                    6: 'ST14207',
                    9: 'ST14195',
                    18: 'ST14205',
                    19: 'ST14186',
                    27: 'ST14212',
                    28: 'ST14211'
                    }
df['strain'] = df['design'].map(lambda x:comp_exp_links[x])

# Pivot it to be like a matrix 
df_pivot = df.pivot(index='model', columns='strain', values='glc_yield_inc')
myColors = ((0.8, 0.0, 0.0, 1.0), (1.0, 1.0, 1.0, 1.0), (0.0, 0.0, 0.8, 1.0))
cmap = LinearSegmentedColormap.from_list('Custom', myColors, 50)
sns.heatmap(df_pivot - 1, cmap=cmap, annot=True, fmt=f'.1f', annot_kws={"size": 9}, center=0,
vmin=-0.5, vmax=0.5)

# Ensure all xticks are displayed
ax = plt.gca()
ax.set_xticks(np.arange(10) + 0.5)  # 10 designs!
ax.set_xticklabels(np.arange(1,11), fontsize=12, rotation=0)  # Set all x-tick labels

# Set the size of the xtick and ytick labels
plt.yticks(fontsize=12)  # Adjust the size as needed for ytick labels
cbar = plt.gcf().axes[-1]  # Get the colorbar axis
cbar.tick_params(labelsize=12)  # Set the tick label size for the colorbar
current_xlabel = plt.gca().get_xlabel()
current_ylabel = plt.gca().get_ylabel()
# plt.xlabel(current_xlabel, fontsize=16)  # Change the font size as needed
plt.xlabel('Design', fontsize=16)  # Change the font size as needed
plt.ylabel(current_ylabel, fontsize=16)  # Change the font size as needed
plt.tight_layout()
plt.savefig('./../../output/figures/Figure_2A.png', bbox_inches='tight')
plt.close()
