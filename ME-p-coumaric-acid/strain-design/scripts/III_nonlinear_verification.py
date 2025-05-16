"""
This script verifies the top 10 robust designs in nonlinear simulations in each of the 9 kinetic models.
"""
import numpy as np
import pandas as pd
import glob

from pytfa.io.json import load_json_model
from pytfa.optim.constraints import *

from skimpy.io.yaml import load_yaml_model

from skimpy.analysis.oracle.load_pytfa_solution import load_fluxes, load_concentrations, load_equilibrium_constants
from skimpy.core.parameters import load_parameter_population
from skimpy.core.solution import ODESolutionPopulation
from skimpy.utils.tabdict import TabDict
from skimpy.utils.general import get_stoichiometry
from skimpy.utils.namespace import NET, QSSA
from skimpy.simulations.reactor import make_batch_reactor

import matplotlib.pyplot as plt
import os
import time


# Do we want to plot all the solutions as well?
PLOT = True

# IX_START varies from 0 - 8 for the 9 models. Here, each time we can go model by model or if we do IX_START = 0
# and IX_END = 9, we will do all the models. but this is computationally very time consuming.
IX_START = 3
IX_END = IX_START + 1

# Reactor parameters
TOTAL_TIME = 48 # hours
N_STEPS = 1000
INIT_BIOMASS = 0.024 # g - corresponds to a starting OD600 of 0.2

# Cellular parameters
CONCENTRATION_SCALING = 1e9 # 1 mol to 1 mmol
TIME_SCALING = 1 # 1hr
DENSITY = 1200 # g/L
GDW_GWW_RATIO = 0.3 # Assumes 70% Water
CELL_VOLUME = 42 # mum3
GDW_PER_CELL = GDW_GWW_RATIO * DENSITY * CELL_VOLUME * 1e-15 # 1e-15 L per mum3
flux_scaling_factor = 1e-3 * (GDW_GWW_RATIO * DENSITY) * CONCENTRATION_SCALING / TIME_SCALING

# Paths
path_to_tmodel = './../../models/constraint_based_model_scaffold.json'
path_to_kmodel = './../../models/kinetic_model_scaffold.yaml'
path_to_tfa_sample = './../../data/steady_states.csv'
path_to_params = './../../kinetic-modelling/parameters/parameters_top_9_models.hdf5'
path_to_designs = './../output/designs/nra_sol_all_designs_all_models.csv'

# Note please change this as you wish - maybe you want to do many studies and want a better naming structure
base_folder_for_output = './../output/nonlinear-verification/'

# Load models, samples and parameters
tmodel = load_json_model(path_to_tmodel)
kmodel = load_yaml_model(path_to_kmodel)
tfa_sample = pd.read_csv(path_to_tfa_sample, index_col=0).iloc[3191] # steady state closest to mean

# Get kinetic parameter sets
parameter_population = load_parameter_population(path_to_params)
kinetic_models = list(parameter_population._index.keys())

# Load all designs with their NRA suggested fold changes and predicted increase in pca yield
df_designs = pd.read_csv(path_to_designs, header=0, index_col=0)

# Get the top 10 designs
top_design_indices = df_designs.groupby('design_ix').mean().sort_values(by='nra_sol', ascending=False).index[:10]
top_designs = df_designs[df_designs['design_ix'].isin(top_design_indices)]

# Prepare the kinetic models 
kmodel.prepare(mca=False)

# Reactor initialization
reactor = make_batch_reactor('single_species.yaml')
reactor.compile_ode(add_dilution=False, ncpu=4)
reactor_volume = reactor.models.strain_1.parameters.strain_1_volume_e.value

concentrations_to_plot = {'glc__D_e': 180 / CONCENTRATION_SCALING,
                    'pca_e': 164 / CONCENTRATION_SCALING,
                    'biomass_strain_1': GDW_PER_CELL }


for this_model in kinetic_models[IX_START: IX_END]:
    print(this_model)
    tfa_ix, kin_ix = this_model.split(',')

    folder_for_output = base_folder_for_output + '{}/'.format(this_model)
    if not os.path.exists(folder_for_output):
        os.makedirs(folder_for_output)

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
        reactor.initial_conditions['biomass_strain_1'] = INIT_BIOMASS * reactor_volume / 1e15 / GDW_PER_CELL # get the actual number of cells from the reactor volume and GDW per cell

        for met_ in reactor.medium.keys():
            if met_.startswith('_'):
                LC_id = 'LC_' + met_[1:]
            else:
                LC_id = 'LC_' + met_
            reactor.initial_conditions[met_] = np.exp(tfa_sample.loc[LC_id]) * 1e9 # 1e9 is the Concentration scaling for numerical purposes.

        # Make sure the reactor volume parameters reflect those of the kinetic model
        for v in [p for p in kmodel.parameters if 'volume' in p]:
            reactor.models.strain_1.parameters['strain_1_' + v].value = kmodel.parameters[v].value
        reactor.models.strain_1.parameters.strain_1_volume_e.value = reactor_volume

    # Root function for stopping criteria
    def rootfn(t, y, g, user_data):
        c_min = user_data['c_min']
        t_0 = user_data['time_0']
        t_max = user_data['max_time']
        t = time.time()
        if (t - t_0) >= t_max:
            g[0] = 0
            print('Oops - did not converge in time')
        else:
            g[0] = 1

    reset_reactor()

    if hasattr(reactor, 'solver'):
        delattr(reactor, 'solver')

    t0 = time.time()
    user_data = {'time_0': t0,
                 'max_time': 36000,  # 10 hours!! The system is quite stiff and it takes a while to integrate!
                 'c_min': 1e-52, }  # If any concentration reaches this low, exit

    TOUT = np.linspace(0, TOTAL_TIME, N_STEPS)
    sol_ode_wt = reactor.solve_ode(TOUT,
                                   solver_type='cvode',
                                   rtol=1e-9,
                                   atol=1e-12,
                                   max_steps=1e9,
                                   rootfn=rootfn,
                                   nr_rootfns=1,
                                   user_data=user_data,
                                   )

    """    
    Test each unique design in a reactor setting for this kinetic model    
    """
    # Get all the enzyme fold changes specific to this model across all 10 designs.
    this_model_designs = top_designs[top_designs['model'] == this_model]
    design_sols = []
    indices = []

    for ix_design, this_design in this_model_designs.groupby('design_ix'):
        reset_reactor()
        indices.append(ix_design)

        # Manipulate the enzymes according to the design
        for e, v in zip(this_design['enzyme'], this_design['fold-change']):            
            parameter_name = "vmax_forward_" + e[3:]
            if 'ED_' in e:
                reactor.parameters['strain_1_' + parameter_name].value *= np.exp(-v)
            elif 'EU_' in e:
                reactor.parameters['strain_1_' + parameter_name].value *= np.exp(v)

        print('Solving design {}'.format(ix_design))

        import time
        if hasattr(reactor, 'solver'):
            delattr(reactor, 'solver')
        t0 = time.time()
        user_data = {'time_0': t0,
                     'max_time': 36000,  # 10 hours!
                     'c_min': 1e-52, }  # If any concentration reaches this low, exit

        sol_ode = reactor.solve_ode(np.linspace(0, TOTAL_TIME, N_STEPS),
                                    solver_type='cvode',
                                    rtol=1e-9,
                                    atol=1e-12,
                                    max_steps=1e9,
                                    rootfn=rootfn,
                                    nr_rootfns=1,
                                    user_data=user_data,
                                    )
        design_sols.append(sol_ode)

        if PLOT == True:
            # Plot this design vs wildtype
            for conc, scaling in concentrations_to_plot.items():
                plt.plot(sol_ode_wt.time, sol_ode_wt.concentrations[conc] * scaling, label='ST10284', color='blue')
                plt.plot(sol_ode.time, sol_ode.concentrations[conc] * scaling, color='orange', label='Design {}'.format(ix_design))
                plt.legend(fontsize=14)
                plt.xlabel('Time (h)', fontsize=14)
                plt.ylabel(conc, fontsize=14)
                plt.xlim([0, TOTAL_TIME])
                plt.savefig(folder_for_output + 'design_{}_{}.png'.format(ix_design, conc), bbox_inches='tight')
                plt.close()

    # Save the ode solutions for the wild type model and all the engineering strains
    solpop = ODESolutionPopulation(design_sols, indices)
    solpop.data.to_csv(folder_for_output + 'sol_designs.csv')
    sol_ode_wt.concentrations.to_csv(folder_for_output + 'sol_wt.csv')