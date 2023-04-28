"""
In this script, we validate the top design for improving the aroG tktA insilico strain in bioreactor simulations
- We use the mean fold changes obtain in the script S.VIII_ii and implement them.
- The top design targets DDPA, TKT1, TKT2, DHQS, SHKK, ANS, and GLUDy
Ouputted data is used for plotting figure S.11
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


'''
Run parameters & paths
'''

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
path_to_tmodel = './../models/tfa_model.json'
path_to_kmodel = './../models/kinetic_model_scaffold.yaml'
path_to_tfa_samples = './../data/tfa_samples.csv'
path_to_params = './../data/kinetic_params_top_10_models.hdf5'
path_to_regulation_data = './../data/allosteric_regulations.csv'

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

# Load set of parameters corresponding to the 10 kinetic models
parameter_population = load_parameter_population(path_to_params)

# Reactor initialization
reactor = make_batch_reactor('single_species.yaml', df_regulation= df_regulations_all)
reactor.compile_ode(add_dilution=False)
reactor_volume = reactor.models.strain_1.parameters.strain_1_volume_e.value

kinetic_models = list(parameter_population._index.keys())

'''
Strain design
1. Go through each kinetic model (combination of a tfa steady state sample and kinetic parameters built around it)
2. Get the control coefficients
3. Generate design alternatives using the steady-state profiles anc control coefficients for a given objective
4. Conduct reactor simulations
5. Store plots and outputs
'''

for this_model in kinetic_models:

    # Get TFA and kinetic model indices
    list_ix = this_model.split(',')
    tfa_ix = int(list_ix[0])

    # Create folder for output
    output_folder = './../output/data/S.VIII-aroGtktA-analysis/reactor-verification/{}/'.format(this_model)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

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
    ''' 
    Reactor simulations 
    '''
    def reset_reactor():

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

    ''' 1. Simulation of the reference strain'''
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


    ''' 3. Paper designs - DDPA'''
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

    ''' 4. Paper designs - DDPA + TKT '''
    reset_reactor()

    # Remove the regulation of DDPA by phenylalanine
    reactor.parameters['strain_1_ki_inhibitor1_DDPA'].value = 1e42
    reactor.parameters['strain_1_k_inhibition_IM_phe_L_c_DDPA'].value = 1e42

    # Upregulate TKT1 and TKT2 by 5 fold
    reactor.parameters['strain_1_vmax_forward_TKT1'].value *= 5
    reactor.parameters['strain_1_vmax_forward_TKT2'].value *= 5

    start = T.time()
    sol_ode_ddpa_tkt = reactor.solve_ode(np.linspace(0, TOTAL_TIME, N_STEPS),
                                    solver_type='cvode',
                                    rtol=1e-9,
                                    atol=1e-9,
                                    max_steps=1e9,
                                    )
    end = T.time()
    print("DDPA + TKT: {}".format(end - start))


    ''' 5. NRA suggested double mutant - DDPA + TKT + shiki'''
    reset_reactor()

    # Remove the regulation of DDPA by phenylalanine
    reactor.parameters['strain_1_ki_inhibitor1_DDPA'].value = 1e42
    reactor.parameters['strain_1_k_inhibition_IM_phe_L_c_DDPA'].value = 1e42

    # Upregulate TKT1 and TKT2 by 5 fold
    reactor.parameters['strain_1_vmax_forward_TKT1'].value *= 4.83
    reactor.parameters['strain_1_vmax_forward_TKT2'].value *= 4.83

    # Upregulate additional enzymes identified by NRA, including DDPA
    reactor.parameters['strain_1_vmax_forward_DDPA'].value *= 5
    reactor.parameters['strain_1_vmax_forward_DHQS'].value *= 5
    reactor.parameters['strain_1_vmax_forward_SHKK'].value *= 5
    reactor.parameters['strain_1_vmax_forward_ANS'].value *= 4.5
    reactor.parameters['strain_1_vmax_forward_GLUDy'].value /= 1.23

    start = T.time()
    sol_ode_ddpa_tkt_shiki = reactor.solve_ode(np.linspace(0, TOTAL_TIME, N_STEPS),
                                    solver_type='cvode',
                                    rtol=1e-9,
                                    atol=1e-9,
                                    max_steps=1e9,
                                    )
    end = T.time()
    print("DDPA + TKT + SHIKI: {}".format(end - start))

    # Store all ode solutions
    sol_odes_all = [sol_ode_wildtype, sol_ode_ddpa, sol_ode_ddpa_tkt, sol_ode_ddpa_tkt_shiki]
    solpop = ODESolutionPopulation(sol_odes_all, np.arange(0, len(sol_odes_all)))
    solpop.data.to_csv(output_folder + 'ode_solutions.csv')

