"""
This is the script to simulate the nonlinear behavior of the 35 models that are calibrated on the behavior of the
reference strain, trpD9923 and the recombinant strain, trpD9923/pJLaroGfbr.
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


''' Run parameters and paths '''
# Cellular parameters
CONCENTRATION_SCALING = 1e9  # 1 mol to 1 mmol
TIME_SCALING = 1  # 1hr
DENSITY = 1105  # g/L
GDW_GWW_RATIO = 0.3  # Assumes 70% Water
flux_scaling_factor = 1e-3 * (GDW_GWW_RATIO * DENSITY) * CONCENTRATION_SCALING / TIME_SCALING

# Ode simulation parameters
TOTAL_TIME = 60
N_STEPS = 1000

# Paths
path_to_tmodel = '../models/tfa_model.json'
path_to_kmodel = '../models/kinetic_model_scaffold.yaml'
path_to_tfa_samples = '../data/steady_state_samples.csv'
path_to_params = '../data/top_35_kinetic_models.hdf5'
path_to_regulation_data = '../data/allosteric_regulations.csv'
output_folder = './../output/data/top-35-models/'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

'''
Model and simulation initiation + preparation
'''

# Load models
tmodel = load_json_model(path_to_tmodel)
kmodel_draft = load_yaml_model(path_to_kmodel)

# Add regulation data to kinetic model
df = pd.read_csv(path_to_regulation_data)
df_regulations_all = df[df['reaction_id'].isin(list(kmodel_draft.reactions.keys()))]
df_regulations_all = df_regulations_all[df_regulations_all['regulator'].isin(list(kmodel_draft.reactants.keys()))]
kmodel = load_enzyme_regulation(kmodel_draft, df_regulations_all)

# Load set of parameters corresponding to the 35 kinetic models
parameter_population = load_parameter_population(path_to_params)

# Reactor initialization
reactor = make_batch_reactor('single_species.yaml', df_regulation= df_regulations_all)
reactor.compile_ode(add_dilution=False)
reactor_volume = reactor.models.strain_1.parameters.strain_1_volume_e.value

kinetic_models = list(parameter_population._index.keys())

def reset_reactor():
    """
    Function to reset the reactor and load the concentrations and parameters before each simulation
    """
    # Parameterize the rector and initialize with incolum and medium data
    reactor.parametrize(parameter_set, 'strain_1')
    reactor.initialize(concentrations, 'strain_1')
    reactor.initial_conditions['biomass_strain_1'] = 0.037 * 0.05 / 0.28e-12

    for met_ in reactor.medium.keys():
        LC_id = 'LC_' + met_
        LC_id = LC_id.replace('_L', '-L')
        LC_id = LC_id.replace('_D', '-D')
        reactor.initial_conditions[met_] = np.exp(tfa_sample.loc[LC_id]) * 1e9

    # Volume parameters for the reactor
    reactor.models.strain_1.parameters.strain_1_volume_e.value = reactor_volume
    reactor.models.strain_1.parameters.strain_1_cell_volume_e.value = 1.0  # 1.0 #(mum**3)
    reactor.models.strain_1.parameters.strain_1_cell_volume_c.value = 1.0  # 1.0 #(mum**3)
    reactor.models.strain_1.parameters.strain_1_cell_volume_p.value = 1.0  # 1.0 #(mum**3)
    reactor.models.strain_1.parameters.strain_1_volume_c.value = 0.9 * 1.0  # (mum**3)
    reactor.models.strain_1.parameters.strain_1_volume_p.value = 0.1 * 1.0  # (mum**3)

sols_wt, sols_d1, sols_d2 = [], [], []
for this_model in kinetic_models:

    # Get TFA and kinetic model indices
    list_ix = this_model.split(',')
    tfa_ix = int(list_ix[0])

    # Load tfa sample and kinetic parameters into kinetic model
    tfa_sample = pd.read_csv(path_to_tfa_samples, header=0, index_col=0).iloc[tfa_ix]
    parameter_set = parameter_population[this_model]
    kmodel.parameters = parameter_set

    # Load steady-state fluxes and concentrations into the scaffold kmodel
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

    ''' 1. Simulation of trpD9923'''
    reset_reactor()
    start = T.time()
    sol_ode_wildtype = reactor.solve_ode(np.linspace(0, TOTAL_TIME, N_STEPS),
                                         solver_type='cvode',
                                         rtol=1e-9,
                                         atol=1e-9,
                                         max_steps=1e9,
                                         )
    end = T.time()
    print("WT: {}".format(end - start))
    sols_wt.append(sol_ode_wildtype)

    ''' 2. Simulation of trpD9923/pJLaroGfbr'''
    reset_reactor()

    # Remove the regulation of DDPA by phenylalanine
    reactor.parameters['strain_1_ki_inhibitor1_DDPA'].value = 1e42
    reactor.parameters['strain_1_k_inhibition_IM_phe_L_c_DDPA'].value = 1e42

    start = T.time()
    sol_ode_ddpa = reactor.solve_ode(np.linspace(0, TOTAL_TIME, N_STEPS),
                                     solver_type='cvode',
                                     rtol=1e-9,
                                     atol=1e-9,
                                     max_steps=1e9,
                                     )
    end = T.time()
    print("DDPA: {}".format(end - start))
    sols_d1.append(sol_ode_ddpa)

    ''' 3. Simulation of trpD9923/pJLaroGfbrtktA'''
    reset_reactor()

    # Remove the regulation of DDPA by phenylalanine
    reactor.parameters['strain_1_ki_inhibitor1_DDPA'].value = 1e42
    reactor.parameters['strain_1_k_inhibition_IM_phe_L_c_DDPA'].value = 1e42

    # Upregulate TKT1 and TKT2 by 10 fold
    reactor.parameters['strain_1_vmax_forward_TKT1'].value *= 10
    reactor.parameters['strain_1_vmax_forward_TKT2'].value *= 10

    start = T.time()
    sol_ode_ddpa_tkt = reactor.solve_ode(np.linspace(0, TOTAL_TIME, N_STEPS),
                                         solver_type='cvode',
                                         rtol=1e-9,
                                         atol=1e-9,
                                         max_steps=1e9,
                                         )
    end = T.time()
    print("DDPA + TKT: {}".format(end - start))
    sols_d2.append(sol_ode_ddpa_tkt)

# Store all ode solutions
solpop = ODESolutionPopulation(sols_wt, kinetic_models)
solpop.data.to_csv(output_folder + 'ode_sols_wt.csv')

solpop = ODESolutionPopulation(sols_d1, kinetic_models)
solpop.data.to_csv(output_folder + 'ode_sols_d1.csv')

solpop = ODESolutionPopulation(sols_d2, kinetic_models)
solpop.data.to_csv(output_folder + 'ode_sols_d2.csv')