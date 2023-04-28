"""
In this script, we apply each pf the top 5 designs to each of the models and get the NRA suggested values.
We also obtain the fold changes for each of the top 5 designs suggested for each kinetic model

NOTE: This script is needed to get Figure 7
"""
import numpy as np
import pandas as pd
import glob

from pytfa.io.json import load_json_model
from pytfa.optim.constraints import *

from skimpy.io.yaml import load_yaml_model
from skimpy.io.regulation import load_enzyme_regulation
from skimpy.analysis.oracle.load_pytfa_solution import load_fluxes, load_concentrations, load_equilibrium_constants
from skimpy.core.parameters import load_parameter_population
from skimpy.utils.tabdict import TabDict
from skimpy.utils.namespace import NET, QSSA
from skimpy.simulations.reactor import make_batch_reactor

from NRAplus.core.model import NRAmodel


DESIGN_INDICES = []

# Design parameters
N_ENZYMES = 3
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

# Load models, samples and parameters
tmodel = load_json_model(path_to_tmodel)
kmodel_draft = load_yaml_model(path_to_kmodel)

# Load the dataframe with regulations and choose those rxns and mets that are in the kinetic model
df = pd.read_csv(path_to_regulation_data)
df_regulations_all = df[df['reaction_id'].isin(list(kmodel_draft.reactions.keys()))]
df_regulations_all = df_regulations_all[df_regulations_all['regulator'].isin(list(kmodel_draft.reactants.keys()))]

# Create kinetic model with regulation
kmodel = load_enzyme_regulation(kmodel_draft, df_regulations_all)

parameter_population = load_parameter_population(path_to_params)
kinetic_models = list(parameter_population._index.keys())

# These are the top 5 designs that were identified
list_of_designs = [['EU_DDPA', 'ED_GLUDy', 'EU_PYK',],# d-1 - circles
                   ['EU_DDPA', 'ED_GLUDy', 'ED_PGI',], # d-2 - dashed line
                   ['EU_DDPA', 'ED_GLUDy', 'ED_GND',],# d-3 - dash dot
                   ['EU_DDPA', 'ED_GLUDy', 'ED_HEX1',],# d-4 - dotted
                   ['EU_DDPA', 'ED_GLUDy', 'EU_ANS',], # d-5 - solid line
                   ]

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
kmodel.compile_mca(mca_type=NET, sim_type=QSSA, parameter_list=parameter_list, ncpu=1)

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

    # Get control coefficients
    flux_control_coeff = kmodel.flux_control_fun(fluxes, concentrations, [parameter_set])
    concentration_control_coeff = kmodel.concentration_control_fun(fluxes, concentrations, [parameter_set])
    ccc = concentration_control_coeff.mean('sample')
    fcc = flux_control_coeff.mean('sample')

    def reset_nmodel(tmodel, tfa_sample, fcc, ccc):
        """ function for repeated calls to create nmodel
        """

        # Remove bounds on anthranilate concentrations in the tmodel first
        for this_LC in tmodel.log_concentration:
            this_LC.variable.ub = 100
            this_LC.variable.lb = -100

        # Build NRAmodel based on the mean control coeffcients
        nmodel = NRAmodel(tmodel, tfa_sample, fcc, ccc, type=NET)

        # Set Objective maximise the relative change of anthranilate secretion
        nmodel.objective = nmodel.variables.LFR_ANS - nmodel.variables.LFR_GLCtex

        # Ensure that TKT is expressed together
        expr = nmodel.variables.EU_TKT1 - nmodel.variables.EU_TKT2
        nmodel.add_constraint(kind=ModelConstraint, id_='TKT_ratio_up_', hook=nmodel, expr=expr, lb=0, ub=0)
        expr = nmodel.variables.ED_TKT1 - nmodel.variables.ED_TKT2
        nmodel.add_constraint(kind=ModelConstraint, id_='TKT_ratio_down_', hook=nmodel, expr=expr, lb=0, ub=0)

        # Design constraints
        nmodel.max_enzyme_modifications = N_ENZYMES
        nmodel.max_fold_enzyme_change = np.log(MAX_FOLD_ENZ_ACTIVITY)
        nmodel.max_fold_concentration_change = np.log(MAX_FOLD_CONC_CHANGE)
        for var in nmodel.enzyme_down_regulation_variable:  # Allow for knockouts
            var.ub = 100
        nmodel.variables.EU_LMPD_biomass_c_1_420.ub = 1e-3  # Don't allow enzymatic changes to the lumped reaction
        nmodel.variables.ED_LMPD_biomass_c_1_420.ub = 1e-3

        nmodel.variables.LFR_LMPD_biomass_c_1_420.lb = -0.223  # Corresponding to 80% of original biomass

        # Numerical configurations
        nmodel.solver.configuration.tolerances.feasibility = 1e-9
        nmodel.solver.configuration.tolerances.optimality = 1e-9
        nmodel.solver.configuration.tolerances.integrality = 1e-9
        nmodel.solver.problem.parameters.read.scale = -1
        nmodel.solver.problem.parameters.emphasis.numerical = 1

        nmodel.solver.configuration.presolve = True #*************** VVVIMP ******************

        return nmodel

    nmodel = reset_nmodel(tmodel, tfa_sample, fcc, ccc)

    """    
    Test each unique design for this kinetic model using NRA    
    """
    this_nmodel = nmodel.copy()
    for ix, this_design in enumerate(list_of_designs):

        # Enforce current design (enzymes on/off)
        for this_enzyme in this_design:
            this_nmodel.variables[this_enzyme].lb = 1e-3 # Set to a minimum level to ensure it is "ON"
        try:
            sol = this_nmodel.optimize()
            obj = sol.objective_value
        except:
            print("ix {} did not work".format(ix))
            obj = -100


        # Store all values in the dictionary
        for this_enzyme in this_design:
            dict_ = {'model': this_model,
                     'design_ix': ix,
                     'Enzyme': this_enzyme,
                     'Fold-change': sol.raw.loc[this_enzyme],
                     'nra_sol': obj,
                     }
            list_dicts.append(dict_)

        # Reset enzyme regulations
        for this_enzyme in this_design:
            this_nmodel.variables[this_enzyme].lb = 0

df_ = pd.DataFrame(list_dicts)
df_.to_csv('../output/data/{}-fold-conc-change/top_5_designs_fold_changes.csv'.format(MAX_FOLD_CONC_CHANGE))


