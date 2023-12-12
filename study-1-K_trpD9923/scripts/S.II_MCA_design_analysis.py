"""
This code is used to run the experiments presented in Supplementary Note II - MCA design.
For each of the 10 models, we do the following
1. We identify the yield control coefficients (YCCs) for anthranilate wrt glucose
2. We then prune the YCCs for those enzymes that have the same sign of the YCC as they do for their
biomass control coefficient. This ensures that perturbing this enzyme does not negatively impact
growth
3. We then choose the top pruned YCCs and perturb them in bioreactor simulations with 1.2 fold,
2 fold, and 5 fold changes in enzyme activity
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


SIMULATE_IN_REACTOR = True

# Design parameters
MAX_FOLD_ENZ_ACTIVITY = [1.2, 2, 5]
N_ENZ = 3

# Cellular parameters
CONCENTRATION_SCALING = 1e9  # 1 mol to 1 mmol
TIME_SCALING = 1  # 1hr
DENSITY = 1105  # g/L
GDW_GWW_RATIO = 0.3  # Assumes 70% Water
flux_scaling_factor = 1e-3 * (GDW_GWW_RATIO * DENSITY) * CONCENTRATION_SCALING / TIME_SCALING
concentrations_to_plot = { 'anth_e': 136.13 / CONCENTRATION_SCALING,
                          'biomass_strain_1': 0.28e-12 / 0.05}

# Ode simulation parameters
TOTAL_TIME = 65
N_STEPS = 1000

# Paths
path_to_tmodel = './../models/tfa_model.json'
path_to_kmodel = './../models/kinetic_model_scaffold.yaml'
path_to_tfa_samples = './../data/tfa_samples.csv'
path_to_params = './../data/kinetic_params_top_10_models.hdf5'
path_to_regulation_data = './../data/allosteric_regulations.csv'


# Load models
tmodel = load_json_model(path_to_tmodel)
kmodel_draft = load_yaml_model(path_to_kmodel)

# Add regulation data to kinetic model
df = pd.read_csv(path_to_regulation_data)
df_regulations_all = df[df['reaction_id'].isin(list(kmodel_draft.reactions.keys()))]
df_regulations_all = df_regulations_all[df_regulations_all['regulator'].isin(list(kmodel_draft.reactants.keys()))]
kmodel = load_enzyme_regulation(kmodel_draft, df_regulations_all)

# Get the correct parameter sets
parameter_population = load_parameter_population(path_to_params)
kinetic_models = list(parameter_population._index.keys())

# Get the list of parameters for which we want control coefficients
transport_reactions = ['vmax_forward_' + r.id for r in tmodel.reactions
                       if (len(r.compartments) > 1)
                       and not ('í' in r.compartments)
                       and not r.id.startswith('LMPD_')
                       and r.id in kmodel.reactions]
transport_reactions.append('vmax_forward_ATPM')
parameter_list = TabDict([(k, p.symbol) for k, p in kmodel.parameters.items()
                          if p.name.startswith('vmax_forward')
                          and str(p.symbol) not in transport_reactions])

# Prepare the kinetic models for MCA
kmodel.prepare()
kmodel.compile_mca(mca_type=NET, sim_type=QSSA, parameter_list=parameter_list, ncpu=4)

# Reactor initialization
reactor = make_batch_reactor('single_species.yaml', df_regulation= df_regulations_all)
reactor.compile_ode(add_dilution=False)
reactor_volume = reactor.models.strain_1.parameters.strain_1_volume_e.value

# Loop through all kinetic models, get NRA designs and test everything in a bioreactor setup
for this_model in kinetic_models:

    # Get TFA and kinetic model indices
    list_ix = this_model.split(',')
    tfa_ix = int(list_ix[0])
    kin_ix = int(list_ix[1])

    folder_for_output = './../output/mca-study/{}'.format(this_model)
    if not os.path.exists(folder_for_output):
        os.makedirs(folder_for_output)

    # Load tfa sample and kinetic parameters into kinetic model
    tfa_sample = pd.read_csv(path_to_tfa_samples, header=0, index_col=0).iloc[tfa_ix]
    parameter_set = parameter_population['{}'.format(this_model)]
    kmodel.parameters = parameter_set

    # Load steady-state fluxes and concentrations into kmodel
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

    """
    GET TOP CCs
    """
    # Get yield control coefficients
    flux_control_coeff = kmodel.flux_control_fun(fluxes, concentrations, [parameter_set])
    concentration_control_coeff = kmodel.concentration_control_fun(fluxes, concentrations, [parameter_set])
    ccc = concentration_control_coeff.mean('sample')
    fcc = flux_control_coeff.mean('sample')
    yield_cc = fcc.loc['ANS'] - fcc.loc['GLCtex']

    # Only those enzymes which don't kill biomass can be targets
    biomass_cc = fcc.loc['LMPD_biomass_c_1_420']

    # Modified yield control coefficients - YCCs that have the same sign for biomass and anthranilate yield
    enz_to_keep = []
    for this_enz in yield_cc.index:
        sign_ycc = np.sign(yield_cc.loc[this_enz])
        sign_gcc = np.sign(biomass_cc.loc[this_enz])

        if sign_ycc == sign_gcc:
            enz_to_keep.append(this_enz) # If they have the same directional impact, keep it

    yield_cc = yield_cc.loc[enz_to_keep]
    highest_ycc = yield_cc.abs().sort_values().index[-N_ENZ:]

    # get combined ccs
    df_combined = pd.concat([yield_cc, biomass_cc.loc[enz_to_keep]], axis=1)
    df_combined.loc[df_combined.abs().sort_values(by=0, ascending=False).index].to_csv(folder_for_output + '/combined_ccs.csv')

    ''' 
    Reactor simulations 
    '''
    def reset_reactor():
        reactor.parametrize(parameter_set, 'strain_1')
        reactor.initialize(concentrations, 'strain_1')
        reactor.initial_conditions['biomass_strain_1'] = 0.037 * 0.05 / 0.28e-12

        for met_ in reactor.medium.keys():
            LC_id = 'LC_' + met_
            LC_id = LC_id.replace('_L', '-L')
            LC_id = LC_id.replace('_D', '-D')
            reactor.initial_conditions[met_] = np.exp(tfa_sample.loc[LC_id]) * 1e9

        # Not nice but necessary
        reactor.models.strain_1.parameters.strain_1_volume_e.value = reactor_volume
        reactor.models.strain_1.parameters.strain_1_cell_volume_e.value = 1.0  # 1.0 #(mum**3) look up cell volume bionumbers
        reactor.models.strain_1.parameters.strain_1_cell_volume_c.value = 1.0  # 1.0 #(mum**3)
        reactor.models.strain_1.parameters.strain_1_cell_volume_p.value = 1.0  # 1.0 #(mum**3)
        reactor.models.strain_1.parameters.strain_1_volume_c.value = 0.9 * 1.0  # (mum**3)
        reactor.models.strain_1.parameters.strain_1_volume_p.value = 0.1 * 1.0  # (mum**3)

    ''' 1. Wildtype '''
    reset_reactor()
    start = T.time()
    sol_ode_wildtype = reactor.solve_ode(np.linspace(0, TOTAL_TIME, N_STEPS),
                                         solver_type='cvode',
                                         rtol=1e-9,
                                         atol=1e-9,
                                         max_steps=1e9,
                                         )
    end = T.time()
    print("WT: {}".format(end-start))

    ''' 2. MCA based designs '''
    sols_ode_fcc = []
    for this_enz_activity in MAX_FOLD_ENZ_ACTIVITY:
        reset_reactor()

        # Manipulate the top fccs proportionally (proportional to their magnitude)
        for this_ycc in highest_ycc:
            reactor.parameters['strain_1_'+this_ycc].value *= \
                this_enz_activity**(np.sign(yield_cc.loc[this_ycc]))

        start = T.time()
        sol_ode_fcc = reactor.solve_ode(np.linspace(0, TOTAL_TIME, N_STEPS),
                                        solver_type='cvode',
                                rtol=1e-9,
                                atol=1e-9,
                                max_steps=1e9,
                                )
        end = T.time()
        print("MCA: {}".format(end-start))
        sols_ode_fcc.append(sol_ode_fcc)

    sol_odes_all = [sol_ode_wildtype]
    sol_odes_all.extend(sols_ode_fcc)
    solpop = ODESolutionPopulation(sol_odes_all, np.arange(0, len(sol_odes_all)))
    solpop.data.to_csv(folder_for_output + '/ode_solutions.csv')
    del solpop

    # Print outputs to file
    yield_cc.loc[highest_ycc].to_csv(folder_for_output + '/top_yccs.csv')
    biomass_cc.loc[highest_ycc].to_csv(folder_for_output + '/top_yccs_growth.csv')
