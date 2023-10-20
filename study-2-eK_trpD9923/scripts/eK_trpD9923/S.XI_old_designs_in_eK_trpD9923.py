'''
In this file we implement the designs generated using K_trpD9923 in eK_trpD9923

We take the top 4 designs from the old models, get their NRA suggested values for each enhanced model under the renewed constraints:
--> allowable fold change in enz act of 10
--> allowable fold change in concs of 2.5
'''
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


# Design parameters
N_ENZYMES = 3
MAX_FOLD_ENZ_ACTIVITY = 10
MAX_FOLD_CONC_CHANGE = 2.5
TOTAL_TIME = 80
N_STEPS = 1000

# Cellular parameters
CONCENTRATION_SCALING = 1e9 # 1 mol to 1 mmol
TIME_SCALING = 1 # 1hr
DENSITY = 1105 # g/L
GDW_GWW_RATIO = 0.3 # Assumes 70% Water
flux_scaling_factor = 1e-3 * (GDW_GWW_RATIO * DENSITY) * CONCENTRATION_SCALING / TIME_SCALING

# Paths
path_to_tmodel = '../../models/tfa_model.json'
path_to_kmodel = '../../models/kinetic_model_scaffold.yaml'
path_to_tfa_samples = '../../data/steady_state_samples.csv'
path_to_params = '../../data/enhanced_kinetic_models.hdf5'
path_to_regulation_data = '../../data/allosteric_regulations.csv'

folder_for_output = './../../output/data/old-designs-in-eK_trpD9923'
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

# Load the top 5 designs obtained using K_trpD9923
list_of_designs = [['EU_DDPA', 'ED_GLUDy', 'EU_PYK',],# d-1 - circles
                   ['EU_DDPA', 'ED_GLUDy', 'ED_PGI',], # d-2 - dashed line
                   ['EU_DDPA', 'ED_GLUDy', 'ED_GND',],# d-3 - dash dot
                   ['EU_DDPA', 'ED_GLUDy', 'ED_HEX1',],# d-4 - dotted
                   ['EU_DDPA', 'ED_GLUDy', 'EU_ANS',], # d-5 - solid line
                   ]

# Get the list of parameters for which we want control coefficients - all vmaxs except transports, ATPM, and biomass
transport_reactions = ['vmax_forward_' + r.id for r in tmodel.reactions
                       if (len(r.compartments) > 1)
                       and not ('í' in r.compartments)
                       and not r.id.startswith('LMPD_')
                       and r.id in kmodel.reactions]
transport_reactions.append('vmax_forward_ATPM')
parameter_list = TabDict([(k, p.symbol) for k, p in kmodel.parameters.items()
                          if p.name.startswith('vmax_forward')
                          and str(p.symbol) not in transport_reactions])

# Prepare the kinetic models for Metabolic Control Analysis (MCA)
kmodel.prepare()
kmodel.compile_mca(mca_type=NET, sim_type=QSSA, parameter_list=parameter_list, ncpu=1)

# Reactor initialization
reactor = make_batch_reactor('single_species.yaml', df_regulation= df_regulations_all)
reactor.compile_ode(add_dilution=False)
reactor_volume = reactor.models.strain_1.parameters.strain_1_volume_e.value

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

def reset_nmodel(tmodel, tfa_sample, fcc, ccc, type=NET):

    nmodel = NRAmodel(tmodel, tfa_sample, fcc, ccc, type=NET)
    # Ensure that TKT is expressed together
    expr = nmodel.variables.EU_TKT1 - nmodel.variables.EU_TKT2
    nmodel.add_constraint(kind=ModelConstraint, id_ = 'TKT_ratio_up', hook=nmodel, expr= expr, lb = 0, ub = 0)
    expr = nmodel.variables.ED_TKT1 - nmodel.variables.ED_TKT2
    nmodel.add_constraint(kind=ModelConstraint, id_ = 'TKT_ratio_down', hook=nmodel, expr= expr, lb = 0, ub = 0)

    # Design constraints
    nmodel.objective = nmodel.variables.LFR_ANS - nmodel.variables.LFR_GLCtex  # anthranilate yield change
    nmodel.max_enzyme_modifications = N_ENZYMES
    nmodel.max_fold_enzyme_change = np.log(MAX_FOLD_ENZ_ACTIVITY)
    nmodel.max_fold_concentration_change = np.log(MAX_FOLD_CONC_CHANGE)

    for var in nmodel.enzyme_down_regulation_variable: # Allow for knockouts
        var.ub = 100

    # Don't allow enzymatic changes to the lumped reaction
    nmodel.variables.EU_LMPD_biomass_c_1_420.ub = 1e-3
    nmodel.variables.ED_LMPD_biomass_c_1_420.ub = 1e-3

    # Numerical configurations
    nmodel.solver.configuration.tolerances.feasibility = 1e-9
    nmodel.solver.configuration.tolerances.optimality = 1e-9
    nmodel.solver.configuration.tolerances.integrality = 1e-9
    nmodel.solver.problem.parameters.read.scale = -1
    nmodel.solver.problem.parameters.emphasis.numerical = 1
    nmodel.solver.configuration.presolve = True  # *************** VVVIMP ******************

    return nmodel

list_dicts = [] # list of dictionaries that will store the different solutions for plotting
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

    # simulate the wildtype response first - just to make sure
    reset_reactor()
    sol_wt = reactor.solve_ode(np.linspace(0, TOTAL_TIME, N_STEPS),
                                solver_type='cvode',
                                rtol=1e-9,
                                atol=1e-9,
                                max_steps=1e9,
                                )

    # Get flux + concentrationc control coefficients and calculate yield control coefficients for anthranilate wrt glc
    fcc = kmodel.flux_control_fun(fluxes, concentrations, [parameter_set]).mean('sample')
    ccc = kmodel.concentration_control_fun(fluxes, concentrations, [parameter_set]).mean('sample')
    
    # Remove bounds on concentrations in the thermo model - we only constrain fold changes in concentrations
    for this_LC in tmodel.log_concentration:
        this_LC.variable.ub = 100
        this_LC.variable.lb = -100

    # Create the NRA model
    nmodel = reset_nmodel(tmodel, tfa_sample, fcc, ccc, type=NET)

    # Loop through the top 5 designs and implement them in this model
    for ix_design, this_design in enumerate(list_of_designs):
        # Enforce current design (enzymes on/off)
        for this_enzyme in this_design:
            nmodel.variables[this_enzyme].lb = 1e-3 # Set to tolerance level
        try:
            sol = nmodel.optimize()
            obj = sol.objective_value
        except:
            print("ix {} did not work".format(ix_design))
            obj = -100

        # Reset enzyme regulations
        for this_enzyme in this_design:
            nmodel.variables[this_enzyme].lb = 0

        # Run reactor simulations of this design for this model
        reset_reactor()
        
        # Manipulate the enzymes according to the NRA design
        for this_enzyme in this_design:
            parameter_name = "vmax_forward_" + this_enzyme[3:]
            if 'ED_' in this_enzyme:
                reactor.parameters['strain_1_'+parameter_name].value *= np.exp(-sol.raw.loc[this_enzyme])
            elif 'EU_' in this_enzyme:
                reactor.parameters['strain_1_'+parameter_name].value *= np.exp(sol.raw.loc[this_enzyme])

        # if DDPA is in the design, dergulate it --> this is the case for the top designs from K_trpD9923
        if 'EU_DDPA' in this_design:
            reactor.parameters['strain_1_ki_inhibitor1_DDPA'].value = 1e42
            reactor.parameters['strain_1_k_inhibition_IM_phe_L_c_DDPA'].value = 1e42

        start = T.time()
        sol_ode = reactor.solve_ode(np.linspace(0, TOTAL_TIME, N_STEPS),
                                solver_type='cvode',
                                rtol=1e-9,
                                atol=1e-9,
                                max_steps=1e9,
                                )

        dict_ = {'model': this_model,
                'design_ix': ix_design,
                'E1': this_design[0],
                'E2': this_design[1],
                'E3': this_design[2],
                'E1_val': sol.raw.loc[this_design[0]],
                'E2_val': sol.raw.loc[this_design[1]],
                'E3_val': sol.raw.loc[this_design[2]],
                'ode_sol_design': sol_ode,
                'ode_sol_wt': sol_wt,
            }
        list_dicts.append(dict_)

# Plot and store!!
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