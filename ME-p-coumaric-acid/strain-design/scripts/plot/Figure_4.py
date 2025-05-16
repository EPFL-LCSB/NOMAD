'''
Script to take batch fermentation data and plot it - Figure 4 in the main manuscript

https://davidmathlogic.com/colorblind/#%23000000-%23E69F00-%2356B4E9-%23009E73-%23F0E442-%230072B2-%23D55E00-%23CC79A7-%2314dff9-%23ea4ac0
'''
import pandas as pd
import matplotlib.pyplot as plt


od = pd.read_csv('./../../../data/batch_experiments_biomass.csv', index_col=0)
pca = pd.read_csv('./../../../data/batch_experiments_pca.csv', index_col=0)

# Colour/marker scheme for designs that performed well and poorly
good_colors = ['blue', 'red', 'gold', 'cyan', 'gray', 'navy', 'maroon', 'brown']
bad_colors = ['blue', 'red',]

good_markers = ['2', 'x', '*', '3', 7, '|', '1', '+']
bad_markers = ['3', 'x']

good_designs = ['ST14211', 'ST14212', 'ST14210', 'ST14186', 'ST14201', 'ST14205', 'ST14206', 'ST14195']
bad_designs = [ 'ST14204', 'ST14207', ]

def plot_curves(df, ax, label='OD600'):
    # First plot all the ones that have better biomass and better pca - all except 14204 and 14207
    for i, d in enumerate(good_designs):
        df_ = df[df.Strain == d]
        ax.errorbar(df_.index, df_.Mean, yerr=df_.SD, color=good_colors[i], label=d, marker='o', markersize=8, capsize=5, capthick=1)

    # next do the ones that were not good
    for i, d in enumerate(bad_designs):
        df_ = df[df.Strain == d]
        ax.errorbar(df_.index, df_.Mean, yerr=df_.SD, color=bad_colors[i], label=d, marker='x', markersize=10, capsize=5, capthick=1)

    df_ = df[df.Strain == 'ST10284']
    ax.errorbar(df_.index, df_.Mean, yerr=df_.SD, linewidth=3, color='black', ls='--', capsize=5, capthick=1)
    ax.set_ylabel(label, fontsize=20)
    ax.tick_params(axis='y', which='major', labelsize=16)

# Create subplots with shared x-axis
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 14))

plot_curves(od, ax2, label='OD600')
plot_curves(pca, ax1, label=' $\it{p}$-CA titer (mg/L)')

ax2.set_xlabel('Time (hours)', fontsize=20)
ax1.legend(loc="lower center", ncol=5, bbox_to_anchor=(0.5,1.1), fontsize=16)

ax2.tick_params(axis='x', which='major', labelsize=16)
ax1.tick_params(labelbottom=False)  # Remove x-axis labels from top plot
ax2.set_xlim([0, 75])
plt.subplots_adjust(hspace=0)  # Adjust space between plots

plt.tight_layout()
plt.savefig('./../../output/figues/Figure_4.png', bbox_inches='tight')
plt.close()


