"""
This script takes in the 10 chosen kinetic models, and generates engineering strains for all of them

USER OPTIONS:
- Change the MAX_FOLD_CONC_CHANGE for different runs (2, 3, 10, 20)
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

# Design parameters
MAX_FOLD_ENZ_ACTIVITY = 5
MAX_FOLD_CONC_CHANGE = 10
MAX_N_ENZ = 3
SIMULATE_IN_REACTOR = True
PLOT_INDIVIDUAL = False

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
    output_folder = './../output/data/{}-fold-conc-change/{}/'.format(MAX_FOLD_CONC_CHANGE, this_model)
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

    """
    NRA analysis
    """
    # Get flux + concentrationc control coefficients and calculate yield control coefficients for anthranilate wrt glc
    fcc = kmodel.flux_control_fun(fluxes, concentrations, [parameter_set]).mean('sample')
    ccc = kmodel.concentration_control_fun(fluxes, concentrations, [parameter_set]).mean('sample')
    yield_cc = fcc.loc['ANS'] - fcc.loc['GLCtex']
    highest_ycc = yield_cc.abs().sort_values().index[-MAX_N_ENZ:]

    # Remove bounds on concentrations in the thermo model - we only constrain fold changes in concentrations
    for this_LC in tmodel.log_concentration:
        this_LC.variable.ub = 100
        this_LC.variable.lb = -100

    # Create the NRA model
    nmodel = NRAmodel(tmodel, tfa_sample, fcc, ccc, type=NET)

    # Ensure that TKT is expressed together
    expr = nmodel.variables.EU_TKT1 - nmodel.variables.EU_TKT2
    nmodel.add_constraint(kind=ModelConstraint, id_ = 'TKT_ratio_up', hook=nmodel, expr= expr, lb = 0, ub = 0)
    expr = nmodel.variables.ED_TKT1 - nmodel.variables.ED_TKT2
    nmodel.add_constraint(kind=ModelConstraint, id_ = 'TKT_ratio_down', hook=nmodel, expr= expr, lb = 0, ub = 0)

    # Design constraints
    nmodel.objective = nmodel.variables.LFR_ANS - nmodel.variables.LFR_GLCtex  # anthranilate yield change
    nmodel.max_enzyme_modifications = MAX_N_ENZ
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

    sol = nmodel.optimize()
    print(sol.objective_value)

    # Enumerate alternatives - max 100 iterations or upto 95% of the solution
    _, designs, _ = enumerate_solutions(nmodel, max_iter=100, threshold=0.95)

    del nmodel

    ''' 
    Reactor simulations 
    '''
    if SIMULATE_IN_REACTOR:

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

        ''' 2. Yield control coefficient based design '''
        reset_reactor()

        # Manipulate the top yccs
        for this_ycc in highest_ycc:
            reactor.parameters['strain_1_'+this_ycc].value *= \
                MAX_FOLD_ENZ_ACTIVITY**(np.sign(yield_cc.loc[this_ycc]))

        start = T.time()
        sol_ode_ycc = reactor.solve_ode(np.linspace(0, TOTAL_TIME, N_STEPS),
                                solver_type='cvode',
                                rtol=1e-9,
                                atol=1e-9,
                                max_steps=1e9,
                                )
        end = T.time()
        print("MCA: {}".format(end-start))

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

        ''' 4. Paper designs - DDPA + TKT'''
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

        ''' 5. Strain designs'''
        sol_ode_alt = []
        for this_design in designs:
            reset_reactor()

            # Manipulate the enzymes according to the NRA design
            for e,v in pd.Series(this_design).iteritems():
                parameter_name = "vmax_forward_" + e[3:]
                if e == 'solution':
                    continue
                if 'ED_' in e:
                    reactor.parameters['strain_1_'+parameter_name].value *= np.exp(-v)
                elif 'EU_' in e:
                    reactor.parameters['strain_1_'+parameter_name].value *= np.exp(v)

            start = T.time()
            sol_ode_alt.append(reactor.solve_ode(np.linspace(0, TOTAL_TIME, N_STEPS),
                                    solver_type='cvode',
                                    rtol=1e-9,
                                    atol=1e-9,
                                    max_steps=1e9,
                                    ))
            end = T.time()
            print("strain : {}".format(end-start))

        if PLOT_INDIVIDUAL:
            """ 
            Plotting each alternative for each model
            """
            concentrations_to_plot = {'glc_D_e': 1/CONCENTRATION_SCALING ,
                                      'anth_e':  136.13/CONCENTRATION_SCALING,
                                      'biomass_strain_1': 0.28e-12 / 0.05}


            for ix_alt in range(len(sol_ode_alt)):

                for conc, scaling in concentrations_to_plot.items():
                    time = np.linspace(0, 80, 1000)
                    plt.plot(time, sol_ode_wildtype.concentrations[conc] * scaling)
                    plt.plot(time, sol_ode_alt[ix_alt].concentrations[conc] * scaling)
                    plt.plot(time, sol_ode_ycc.concentrations[conc] * scaling)
                    plt.plot(time, sol_ode_ddpa.concentrations[conc] * scaling)
                    plt.plot(time, sol_ode_ddpa_tkt.concentrations[conc] * scaling)
                    plt.legend(['wild-type', 'NOMAD-design', 'MCA-design', 'ddpa', 'ddpa+tkt'])
                    plt.xlabel('time (h)')
                    plt.ylabel(conc)
                    plt.savefig(output_folder+'{}_design_{}.png'.format(conc, ix_alt))
                    plt.close()

        # Store all ode solutions
        sol_odes_all = [sol_ode_wildtype, sol_ode_ycc, sol_ode_ddpa, sol_ode_ddpa_tkt]
        sol_odes_all.extend(sol_ode_alt)
        solpop = ODESolutionPopulation(sol_odes_all, np.arange(0, len(sol_odes_all)))
        solpop.data.to_csv(output_folder + 'ode_solutions.csv')

    # Print designs and top yield control coefficients to file
    df_designs = pd.DataFrame(designs)
    cols = df_designs.columns.to_list()
    cols_ = cols[:MAX_N_ENZ] + cols[MAX_N_ENZ+1:]
    cols_.append(cols[MAX_N_ENZ])
    df_designs[cols_].to_csv(output_folder + 'designs.csv')
    yield_cc.loc[highest_ycc].to_csv(output_folder + 'top_yccs.csv')
