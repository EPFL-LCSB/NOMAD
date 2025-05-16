"""
This code takes in the 9 kinetic models and generates combinations of enzyme modifications that will increase pca yield in glucose.
It uses NRA for this.
"""
import numpy as np
import pandas as pd

from pytfa.io.json import load_json_model
from pytfa.optim.constraints import *

from skimpy.io.yaml import load_yaml_model
from skimpy.analysis.oracle.load_pytfa_solution import load_fluxes, load_concentrations, load_equilibrium_constants
from skimpy.analysis.mca.utils import get_dep_indep_vars_from_basis
from skimpy.core.parameters import load_parameter_population

from skimpy.utils.tabdict import TabDict
from skimpy.utils.general import get_stoichiometry
from skimpy.utils.namespace import NET, QSSA

from NRAplus.core.model import NRAmodel
from NRAplus.analysis.design import enumerate_solutions

from sympy import Matrix
from scipy.sparse import csc_matrix as sparse_matrix

import os


# Design parameters
MAX_FOLD_ENZ_ACTIVITY = 2
MAX_FOLD_CONC_CHANGE = 5
N_ENZ = 3

# Cellular parameters
CONCENTRATION_SCALING = 1e9  # 1 mol to 1 mmol
TIME_SCALING = 1  # 1hr
DENSITY = 1200  # g/L
GDW_GWW_RATIO = 0.3  # Assumes 70% Water
flux_scaling_factor = 1e-3 * (GDW_GWW_RATIO * DENSITY) * CONCENTRATION_SCALING / TIME_SCALING

# Paths
path_to_tmodel = './../../models/constraint_based_model_scaffold.json'
path_to_kmodel = './../../models/kinetic_model_scaffold.yaml'
path_to_tfa_sample = './../../data/steady_states.csv'
path_to_params = './../../kinetic-modelling/parameters/parameters_top_9_models.hdf5'
path_to_conservation_relations = '../../data/ST10284_cons_relations_annotated.csv'

# Load models
tmodel = load_json_model(path_to_tmodel)
kmodel = load_yaml_model(path_to_kmodel)
conservations = pd.read_csv(path_to_conservation_relations, index_col=0)
print('*** Models Loaded ***')

# Get the list of parameters for which we want control coefficients - all except transport reactions.
transport_reactions = ['vmax_forward_' + r.id for r in tmodel.reactions
                       if (len(r.compartments) > 1)
                       and not ('í' in r.compartments)
                       and not r.id.startswith('LMPD_')
                       and r.id in kmodel.reactions]
transport_reactions.append('vmax_forward_ATPM')
parameter_list = TabDict([(k, p.symbol) for k, p in kmodel.parameters.items()
                          if p.name.startswith('vmax_forward')
                          and str(p.symbol) not in transport_reactions])

# Prepare the kinetic models for MCA. We first prepare the model without setting it up for MCA as we will add conservation relations manually
kmodel.prepare(mca=False)

# Transform the conservation realtsion to relfect the dependent and independent weights
# (use set of orthogonal dependents).
L0, pivot = Matrix(conservations[kmodel.reactants].values).rref()
kmodel.conservation_relation = sparse_matrix(L0,dtype=np.float)
dep_ix, indep_ix = get_dep_indep_vars_from_basis(kmodel.conservation_relation)
kmodel.independent_variables_ix = indep_ix
kmodel.dependent_variables_ix = dep_ix
kmodel.reduced_stoichiometry = get_stoichiometry(kmodel,kmodel.reactants)[indep_ix,:]

kmodel.compile_mca(mca_type=NET, sim_type=QSSA, parameter_list=parameter_list, ncpu=10)
print('*** MCA compiled ***')

parameter_population = load_parameter_population(path_to_params)
tfa_sample = pd.read_csv(path_to_tfa_sample, index_col=0).iloc[3191]
kinetic_models = parameter_population._index.keys()

# Loop through all kinetic models, get NRA designs and test everything in a bioreactor setup
for ix_ in kinetic_models:

    # Get TFA and kinetic model indices
    TFA_IX, KIN_IX = ix_.split(',')
    TFA_IX = int(TFA_IX)

    # Load parameter population
    parameter_set = parameter_population[ix_]

    # Create folder for output
    output_folder = './../output/designs/{},{}/'.format(TFA_IX, KIN_IX,)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Load tfa sample and kinetic parameters into kinetic model
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
    NRA analysis using NET control coeffcients
    """
    # Get yield control coefficients with 3x the fluxes because the new kinetic models have 3x the original Vmaxs!!
    flux_control_coeff = kmodel.flux_control_fun(fluxes * 3, concentrations, [parameter_set])
    concentration_control_coeff = kmodel.concentration_control_fun(fluxes * 3, concentrations, [parameter_set])
    ccc = concentration_control_coeff.mean('sample')
    fcc = flux_control_coeff.mean('sample')
    yield_cc = fcc.loc['PCAtex'] - fcc.loc['GLCt1']
    highest_ycc = yield_cc.abs().sort_values().index[-N_ENZ:]

    # The DGo of E4Ptm needs to be relaxed
    tmodel.variables['DGo_E4Ptm'].lb = -2

    # Opening up the log concentrations as we will only apply bounds on the *changes* to concentrations.
    for this_LC in tmodel.log_concentration:
        this_LC.variable.lb = -100
        this_LC.variable.ub = 100

    nmodel = NRAmodel(tmodel, tfa_sample, fcc, ccc, type=NET)

    # PCA yield wrt glucose is the objective. Normally we would use the production fluxe instead of the transport flux but in this case
    # there are two routes for production - from tyr and from phe
    nmodel.objective = nmodel.variables.LFR_PCAtex - nmodel.variables.LFR_GLCt1
    nmodel.variables.LFR_LMPD_s_0450_c_1_256.lb = -0.2

    # Design constraints
    nmodel.max_enzyme_modifications = N_ENZ
    nmodel.max_fold_enzyme_change = np.log(MAX_FOLD_ENZ_ACTIVITY)
    nmodel.max_fold_concentration_change = np.log(MAX_FOLD_CONC_CHANGE)
    nmodel.variables.EU_DDPAm.ub = 0 # Not permitting since ARO3 has already been targetted
    nmodel.variables.ED_DDPAm.ub = 0

    # Ensure that TKT reactions are expressed together because both are encoded by transketolase YBR117C or YPR074C
    expr = nmodel.variables.EU_TKT1 - nmodel.variables.EU_TKT2
    nmodel.add_constraint(kind=ModelConstraint, id_ = 'TKT_ratio_up', hook=nmodel, expr= expr, lb = 0, ub = 0)
    expr = nmodel.variables.ED_TKT1 - nmodel.variables.ED_TKT2
    nmodel.add_constraint(kind=ModelConstraint, id_ = 'TKT_ratio_down', hook=nmodel, expr= expr, lb = 0, ub = 0)

    # Numerical configurations
    nmodel.solver.configuration.tolerances.feasibility = 1e-9
    nmodel.solver.configuration.tolerances.optimality = 1e-9
    nmodel.solver.configuration.tolerances.integrality = 1e-9
    nmodel.solver.problem.parameters.read.scale = -1
    nmodel.solver.problem.parameters.emphasis.numerical = 1
    nmodel.solver.configuration.presolve = True  # *************** VVVIMP ******************

    sol = nmodel.optimize()
    print(sol.objective_value)

    # Enumerate alternatives
    design_sols, designs, n_1 = enumerate_solutions(nmodel, max_iter=100, threshold=0.95)

    del nmodel

    # Print outputs to file
    df_designs = pd.DataFrame(designs)
    cols = df_designs.columns.to_list()
    cols_ = cols[:N_ENZ] + cols[N_ENZ+1:]
    cols_.append(cols[N_ENZ])
    df_designs[cols_].to_csv(output_folder + 'designs.csv')
    yield_cc.loc[highest_ycc].to_csv(output_folder + 'top_yccs.csv')


