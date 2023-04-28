'''
Script to generate figure S.9 in supplementary note S.VI

- Loop through all possible combinations of max allowable enz activity and max permissible concentration changes
- For each model in each combination, store the following
-- Get the NRA predicted inc in anthranilate yield
-- Max anthranilate titer from the simulations
-- Time to max titer from the simulations
-- ratio - max_anth_titer / t_max_titer
-- Current allowable enz activity change
-- Current allowable conc activity change
-- model ID
Finally store all this info in a dataframe

- For each allowable enz fold change
-- Calculate mean predicted inc in anthranilate yield across models for each allowable conc fold change
-- Calculate mean production rate (ratio) across models for each allowable conc fold change
-- Plot these means as a function of allowable conc fold change (figure S.9)

'''


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os


plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'figure.autolayout': True})

# Parameters and paths
MAX_FOLD_ENZ_ACTIVITIES = [2, 3, 5, 10]
MAX_FOLD_CONC_CHANGES = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20] # , 3, 5
folder_to_data = './../../output/data/S.VI-phenotype-perturbation-sensitivity/{}-fold-conc-change-{}-fold-enz-activity/{}'
folder_for_output = './../../output/figures/figure-S.9'
kinetic_models = ['1712,6', '1715,6', '1717,6', '2392,6', '2482,8', '3927,1', '4230,7', '4241,5', '4456,7', '4468,6']

if not os.path.exists(folder_for_output):
    os.makedirs(folder_for_output)

dict_info = [] # Dictionary to store all the relevant information

for c_ in MAX_FOLD_CONC_CHANGES:
    for e_ in MAX_FOLD_ENZ_ACTIVITIES:
        for m_ in kinetic_models:

            # Get the ode solutions and NRA predicted inc in anthranilate yield
            ode_sols = pd.read_csv(folder_to_data.format(c_, e_, m_) + '/ode_solutions.csv', index_col=0, header=0)
            design_info = pd.read_csv(folder_to_data.format(c_, e_, m_) + '/designs.csv', index_col=0, header=0)

            # Get max titer and time to max titer
            ode_sol_design = ode_sols[ode_sols['solution_id'] == 1]
            final_anth = np.max(ode_sol_design['anth_e'])
            final_biomass = ode_sol_design.iloc[-1]['biomass_strain_1'] * 0.28e-12 / 0.05
            ix_final_anth = np.where(ode_sol_design['anth_e'] >= 0.99 * final_anth)[0][0]
            t_final_anth = ode_sol_design.iloc[ix_final_anth]['time']

            # Store all the info in a dictionary
            dict_ = {'model_id': m_,
                     'max_conc': c_,
                     'max_enz_act': e_,
                     'nra_sol': design_info.loc[0]['solution'],
                     'final_anth': final_anth * 136.19 * 1e-9,
                     'final_biomass': final_biomass,
                     'time_final_anth': t_final_anth,
                     }
            dict_info.append(dict_)

df_info = pd.DataFrame(dict_info)
df_info['ratio'] = df_info['final_anth'] / df_info['time_final_anth'] # ratio = production rate

''' PLOTTING '''

# 1. For each allowable enz activity, plot how the mean ratio across all models varies as we vary the phenotype proximity
# from close (2-fold conc change) to far (20-fold conc change)
colors = ['blue', 'black', 'red', 'orange']
j=0
for c_, df_ in df_info.groupby('max_enz_act'):

    # Get mean nra sols and ratios for each allowable concentration across models
    df_mean = df_.groupby('max_conc').mean()
    list_concs = list(df_mean.index)

    plt.plot(list_concs, df_mean['ratio'], label= str(c_) + '-fold', color=colors[j])
    j+=1

plt.legend()
plt.xticks(ticks=np.arange(2,20, 2))
plt.xlabel('Allowable fold change in concentrations')
plt.ylabel('Production rate (g/hr)')
plt.savefig(folder_for_output + '/production_rate.png')
plt.close()

# 2. Plot mean anthranilate titers across all models varies as we vary the phenotype proximity, for each allowable enz
# activity change
j=0
for c_, df_ in df_info.groupby('max_enz_act'):

    # Get mean nra sols and ratios for each allowable concentration across models
    df_mean = df_.groupby('max_conc').mean()
    list_concs = list(df_mean.index)

    plt.plot(list_concs, df_mean['final_anth'], label= str(c_) + '-fold', color=colors[j])
    j+=1

plt.legend()
plt.xticks(ticks=np.arange(2,20, 2))
plt.xlabel('Allowable fold change in concentrations')
plt.ylabel('Anthranilate titer (g/L)')
plt.savefig(folder_for_output + '/anthranilate_titers.png')
plt.close()

# 3. Plot mean NRA predicted inc. across all models varies as we vary the phenotype proximity, for each allowable enz
# activity change
j=0
for c_, df_ in df_info.groupby('max_enz_act'):

    # Get mean nra sols and ratios for each allowable concentration across models
    df_mean = df_.groupby('max_conc').mean()
    list_concs = list(df_mean.index)

    plt.plot(list_concs, 100 * (np.exp(df_mean['nra_sol']) - 1), label= str(c_) + '-fold', color=colors[j])
    j+=1

plt.legend()
plt.xticks(ticks=np.arange(2,20, 2))
plt.xlabel('Allowable fold change in concentrations')
plt.ylabel('NRA predicted inc. in anthranilate (%)')
plt.savefig(folder_for_output + '/nra_predictions.png')
plt.close()

# 4. Plot mean growth as we vary the phenotype proximity, for each allowable enz activity change
j=0
for c_, df_ in df_info.groupby('max_enz_act'):

    # Get mean nra sols and ratios for each allowable concentration across models
    df_mean = df_.groupby('max_conc').mean()
    list_concs = list(df_mean.index)

    plt.plot(list_concs, df_mean['final_biomass'], label= str(c_) + '-fold', color=colors[j])
    j+=1

plt.legend()
plt.xticks(ticks=np.arange(2,20, 2))
plt.xlabel('Allowable fold change in concentrations')
plt.ylabel('Biomass titer (g/L)')
plt.savefig(folder_for_output + '/biomass_titers.png')
plt.close()