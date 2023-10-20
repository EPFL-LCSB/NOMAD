"""
We found that DDPA and GLUDy appeared in each of our top 5 designs. Hence, we decided to test the performance of just
this double mutant as suggested by one of the reviewers. The data is used to generate figure S.10 in supplementary note
S.VII
"""
import numpy as np
import pandas as pd

from pytfa.io.json import load_json_model
from pytfa.optim.constraints import *

from skimpy.io.yaml import load_yaml_model
from skimpy.io.regulation import load_enzyme_regulation
from skimpy.analysis.oracle.load_pytfa_solution import load_fluxes, load_concentrations, load_equilibrium_constants
from skimpy.core.parameters import load_parameter_population
from skimpy.core.solution import ODESolutionPopulation
from skimpy.utils.tabdict import TabDict
from skimpy.utils.namespace import NET, QSSA
from skimpy.simulations.reactor import make_batch_reactor

from NRAplus.core.model import NRAmodel
from NRAplus.analysis.design import enumerate_solutions

import os
import glob
import time as T
import matplotlib.pyplot as plt
from sys import argv


# Reactor parameters
TOTAL_TIME = 60
N_STEPS = 1000
MAX_FOLD_ENZ_ACTIVITY = 5
MAX_FOLD_CONC_CHANGE = 3

# Cellular parameters
CONCENTRATION_SCALING = 1e9 # 1 mol to 1 mmol
TIME_SCALING = 1 # 1hr
DENSITY = 1105 # g/L
GDW_GWW_RATIO = 0.3 # Assumes 70% Water
flux_scaling_factor = 1e-3 * (GDW_GWW_RATIO * DENSITY) * CONCENTRATION_SCALING / TIME_SCALING

# Paths
path_to_tmodel = './../models/tfa_model.json'
path_to_kmodel = './../models/kinetic_model_scaffold.yaml'
path_to_tfa_samples = './../data/tfa_samples.csv'
path_to_params = './../data/kinetic_params_top_10_models.hdf5'
path_to_regulation_data = './../data/allosteric_regulations.csv'
folder_for_output = './../output/data/S.VII-double-mutant'
folder_for_figures = './../output/figures/figure-S.10/'

if not os.path.exists(folder_for_output):
    os.makedirs(folder_for_output)

if not os.path.exists(folder_for_figures):
    os.makedirs(folder_for_figures)

# Load models, samples and parameters
tmodel = load_json_model(path_to_tmodel)
kmodel_draft = load_yaml_model(path_to_kmodel)

# Load the dataframe with regulations and choose those rxns and mets that are in the kinetic model
df = pd.read_csv(path_to_regulation_data)
df_regulations_all = df[df['reaction_id'].isin(list(kmodel_draft.reactions.keys()))]
df_regulations_all = df_regulations_all[df_regulations_all['regulator'].isin(list(kmodel_draft.reactants.keys()))]

# Create kinetic model with regulation
kmodel = load_enzyme_regulation(kmodel_draft, df_regulations_all)

# Get the correct parameter sets
parameter_population = load_parameter_population(path_to_params)
kinetic_models = list(parameter_population._index.keys())

# Get the parameters that we want to do MCA for
transport_reactions = ['vmax_forward_'+r.id for r in tmodel.reactions
                       if (len(r.compartments) > 1)
                       and not ('í' in r.compartments)
                       and not r.id.startswith('LMPD_')
                       and r.id in kmodel.reactions]
transport_reactions.append('vmax_forward_ATPM')

parameter_list = TabDict([(k, p.symbol) for k, p in kmodel.parameters.items()
                         if p.name.startswith('vmax_forward')
                         and str(p.symbol) not in transport_reactions])

# Prepare the kinetic model for analysis
kmodel.prepare()

# Reactor initialization
reactor = make_batch_reactor('single_species.yaml', df_regulation= df_regulations_all)
reactor.compile_ode(add_dilution=False)
reactor_volume = reactor.models.strain_1.parameters.strain_1_volume_e.value

list_dicts = []
for this_model in kinetic_models:
    print(this_model)
    tfa_ix, _ = this_model.split(',')

    tfa_sample = pd.read_csv(path_to_tfa_samples, header=0, index_col=0).iloc[int(tfa_ix)]
    parameter_set = parameter_population['{}'.format(this_model)]
    kmodel.parameters = parameter_set

    # Load fluxes and concentrations
    fluxes = load_fluxes(tfa_sample, tmodel, kmodel,
                         density=DENSITY,
                         ratio_gdw_gww=GDW_GWW_RATIO,
                         concentration_scaling=CONCENTRATION_SCALING,
                         time_scaling=TIME_SCALING)

    concentrations = load_concentrations(tfa_sample, tmodel, kmodel,
                                         concentration_scaling=CONCENTRATION_SCALING)

    # Fetch equilibrium constants
    load_equilibrium_constants(tfa_sample, tmodel, kmodel,
                               concentration_scaling=CONCENTRATION_SCALING,
                               in_place=True)

    def reset_reactor():
        reactor.parametrize(parameter_set, 'strain_1')
        reactor.initialize(concentrations, 'strain_1')
        reactor.initial_conditions['biomass_strain_1'] = 0.037 * 0.05 / 0.28e-12

        for met_ in reactor.medium.keys():
            LC_id = 'LC_' + met_
            LC_id = LC_id.replace('_L', '-L')
            LC_id = LC_id.replace('_D', '-D')
            reactor.initial_conditions[met_] = np.exp(tfa_sample.loc[LC_id]) * 1e9

        # Volume settings
        reactor.models.strain_1.parameters.strain_1_volume_e.value = reactor_volume
        reactor.models.strain_1.parameters.strain_1_cell_volume_e.value = 1.0  # 1.0 #(mum**3) look up cell volume bionumbers
        reactor.models.strain_1.parameters.strain_1_cell_volume_c.value = 1.0  # 1.0 #(mum**3)
        reactor.models.strain_1.parameters.strain_1_cell_volume_p.value = 1.0  # 1.0 #(mum**3)
        reactor.models.strain_1.parameters.strain_1_volume_c.value = 0.9 * 1.0  # (mum**3)
        reactor.models.strain_1.parameters.strain_1_volume_p.value = 0.1 * 1.0  # (mum**3)

    reset_reactor()
    sol_ode_wt = reactor.solve_ode(np.linspace(0, TOTAL_TIME, N_STEPS),
                                    solver_type='cvode',
                                    rtol=1e-9,
                                    atol=1e-9,
                                    max_steps=1e9,
                                    )
    """    
    For each model, test with DDPA up regulatd by 4 fold and GLUDy downregulated by 1.2 fold    
    """

    reset_reactor()
    # Manipulate GLUDy and DDPA by their mean values across the top 4 designs
    reactor.parameters['strain_1_vmax_forward_DDPA'].value *= 4 # actual value is 3.9675 - rounded to 4
    reactor.parameters['strain_1_vmax_forward_GLUDy'].value *= 1/1.2 # actual mean is 1.2

    sol_ode = reactor.solve_ode(np.linspace(0, TOTAL_TIME, N_STEPS),
                                solver_type='cvode',
                                rtol=1e-9,
                                atol=1e-9,
                                max_steps=1e9,
                                )

    # Store all values in the dictionary
    dict_ = {'model': this_model,
             'ode_sol_design': sol_ode,
             'ode_sol_wt': sol_ode_wt,
             }
    list_dicts.append(dict_)


df_ = pd.DataFrame(list_dicts)
concentrations_to_plot = {'glc_D_e': 1 / CONCENTRATION_SCALING,
                          'anth_e': 136.13 / CONCENTRATION_SCALING,
                          'biomass_strain_1': 0.28e-12 / 0.05}


# Plotting - NOTE this is now done separately by the script in the plot folder (figure_S.10)
import matplotlib.pyplot as plt

df_sols_design = []
df_sols_wt = []
for ix in range(10):
    df_design = df_.iloc[ix]['ode_sol_design'].concentrations
    df_wt = df_.iloc[ix]['ode_sol_wt'].concentrations
    df_design['time'] = df_.iloc[ix]['ode_sol_design'].time
    df_wt['time'] = df_.iloc[ix]['ode_sol_wt'].time
    df_sols_design.append(df_design)
    df_sols_wt.append(df_wt)

df_sols_wt = pd.concat(df_sols_wt)
df_sols_design = pd.concat(df_sols_design)
#
# df_wt_mean = df_sols_wt.groupby('time').quantile(0.5)
# df_wt_upper = df_sols_wt.groupby('time').quantile(0.75)
# df_wt_lower = df_sols_wt.groupby('time').quantile(0.25)
#
# df_designs_mean = df_sols_design.groupby('time').quantile(0.5)
# df_designs_upper = df_sols_design.groupby('time').quantile(0.75)
# df_designs_lower = df_sols_design.groupby('time').quantile(0.25)

# for conc, scaling in concentrations_to_plot.items():
#     time = df_wt_mean.index
#     # Plot wildtype
#     plt.plot(time, df_wt_mean[conc] * scaling, color='orange', label='wt')
#     plt.fill_between(time, df_wt_lower[conc] * scaling, df_wt_upper[conc] * scaling, facecolor='orange',
#                      interpolate=True,
#                      alpha=0.1)
#
#     # Plot design alternatives
#     plt.plot(time, df_designs_mean[conc] * scaling, color='blue', label='NOMAD - DDPA + GLUDy')
#     plt.fill_between(time, df_designs_lower[conc] * scaling, df_designs_upper[conc] * scaling,
#                      facecolor='blue',
#                      interpolate=True,
#                      alpha=0.1)
#
#     plt.legend()
#     plt.xlabel('time (h)')
#     plt.ylabel(conc.replace('strain_1', '') + ' (g/L)')
#     plt.savefig(folder_for_figures + '/ddpa_gludy_only_{}.png'.format(conc))
#     plt.close()

df_sols_wt.to_csv(folder_for_output + '/ode_sols_wt.csv')
df_sols_design.to_csv(folder_for_output + '/ode_sols_ddpa_gludy.csv')
