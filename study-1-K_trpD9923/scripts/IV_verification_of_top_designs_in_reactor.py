"""
In this script, we verify the top 5 designs in a reactor setup using their model specific NRA suggested fold changes in enzyme activity
"""
import numpy as np
import pandas as pd

from pytfa.io.json import load_json_model
from pytfa.optim.constraints import *

from skimpy.io.yaml import load_yaml_model
from skimpy.io.regulation import load_enzyme_regulation
from skimpy.analysis.oracle.load_pytfa_solution import load_fluxes, load_concentrations, load_equilibrium_constants
from skimpy.core.parameters import load_parameter_population
from skimpy.utils.tabdict import TabDict
from skimpy.simulations.reactor import make_batch_reactor

import os


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
folder_for_output = './../output/verification-top-5-designs-reactor'

if not os.path.exists(folder_for_output):
    os.makedirs(folder_for_output)

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

# Load all NRA designs for all models
list_of_designs = [['EU_DDPA', 'ED_GLUDy', 'EU_PYK',],# d-1 - circles
                   ['EU_DDPA', 'ED_GLUDy', 'ED_PGI',], # d-2 - dashed line
                   ['EU_DDPA', 'ED_GLUDy', 'ED_GND',],# d-3 - dash dot
                   ['EU_DDPA', 'ED_GLUDy', 'ED_HEX1',],# d-4 - dotted
                   ['EU_DDPA', 'ED_GLUDy', 'EU_ANS',], # d-5 - solid line
                   ]

# Get the requisite NRA suggested fold changes for each design for each model
df_enz_act = pd.read_csv('./../output/data/3-fold-conc-change/top_5_designs_fold_changes.csv', header=0, index_col=0)

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
    Test each unique design in a reactor setting for this kinetic model using NRA    
    """

    for ix_design, this_design in enumerate(list_of_designs):

        reset_reactor()

        # Deregulate DDPA since all designs have it as a target
        reactor.parameters['strain_1_ki_inhibitor1_DDPA'].value = 1e42
        reactor.parameters['strain_1_k_inhibition_IM_phe_L_c_DDPA'].value = 1e42

        # Manipulate the enzymes according to the design
        for e in this_design:

            # Get NRA suggested fold change for this enz for this model for this design
            fold_change_row = df_enz_act[(df_enz_act['Enzyme'] == e) &
                                         (df_enz_act['design_ix'] == ix_design) &
                                         (df_enz_act['model'] == this_model)]
            enz_fold_change = fold_change_row['Fold-change'].values[0]

            parameter_name = "vmax_forward_" + e[3:]
            if 'ED_' in e:
                reactor.parameters['strain_1_' + parameter_name].value *= np.exp(-enz_fold_change)
            elif 'EU_' in e:
                reactor.parameters['strain_1_' + parameter_name].value *= np.exp(enz_fold_change)

        print('Solving design {}'.format(ix_design))
        sol_ode = reactor.solve_ode(np.linspace(0, TOTAL_TIME, N_STEPS),
                                    solver_type='cvode',
                                    rtol=1e-9,
                                    atol=1e-9,
                                    max_steps=1e9,
                                    )

        # Store all values in the dictionary
        dict_ = {'model': this_model,
                 'design_ix': ix_design,
                 'E1': this_design[0],
                 'E2': this_design[1],
                 'E3': this_design[2],
                 'ode_sol_design': sol_ode,
                 'ode_sol_wt': sol_ode_wt,
                 }
        list_dicts.append(dict_)


df_ = pd.DataFrame(list_dicts)
concentrations_to_plot = {'glc_D_e': 1 / CONCENTRATION_SCALING,
                          'anth_e': 136.13 / CONCENTRATION_SCALING,
                          'biomass_strain_1': 0.28e-12 / 0.05}

import matplotlib.pyplot as plt

# Group the solutions by design and plot them (wt vs design)
for name, group in df_.groupby('design_ix'):
    df_sols_designs = []
    df_sols_wt = []

    for index, this_row in group.iterrows():
        this_ode_sol_design = this_row['ode_sol_design'].concentrations
        this_ode_sol_wt = this_row['ode_sol_wt'].concentrations
        this_ode_sol_wt['time'] = this_row['ode_sol_wt'].time
        this_ode_sol_design['time'] = this_row['ode_sol_design'].time
        df_sols_designs.append(this_ode_sol_design)
        df_sols_wt.append(this_ode_sol_wt)

    df_sols_designs = pd.concat(df_sols_designs)
    df_sols_wt = pd.concat(df_sols_wt)

    df_sols_designs.to_csv(folder_for_output + '/design_{}_ode_sols.csv'.format(name))

    df_wt_mean = df_sols_wt.groupby('time').quantile(0.5)
    df_wt_upper = df_sols_wt.groupby('time').quantile(0.75)
    df_wt_lower = df_sols_wt.groupby('time').quantile(0.25)

    df_designs_mean = df_sols_designs.groupby('time').quantile(0.5)
    df_designs_upper = df_sols_designs.groupby('time').quantile(0.75)
    df_designs_lower = df_sols_designs.groupby('time').quantile(0.25)

    for conc, scaling in concentrations_to_plot.items():

        time = df_wt_mean.index
        # Plot wildtype
        plt.plot(time, df_wt_mean[conc] * scaling, color='orange', label='wt')
        plt.fill_between(time, df_wt_lower[conc] * scaling, df_wt_upper[conc] * scaling, facecolor='orange',
                         interpolate=True,
                         alpha=0.1)

        # Plot design alternatives
        plt.plot(time, df_designs_mean[conc] * scaling, color='blue', label='nra')
        plt.fill_between(time, df_designs_lower[conc] * scaling, df_designs_upper[conc] * scaling,
                         facecolor='blue',
                         interpolate=True,
                         alpha=0.1)

        plt.legend()
        plt.xlabel('time (h)')
        plt.ylabel(conc.replace('strain_1', '') + ' (g/L)')
        plt.savefig(folder_for_output + '/design_{}_combined_{}.png'.format(name, conc))
        plt.close()

df_sols_wt.to_csv(folder_for_output + '/ode_sols_wt.csv')