"""
In this script, we conduct design robustness analysis
- we take in the 9 kinetic models, and all NRA designs across all models
- we then prune for unique designs (unique by membership, not value of enzyme fold change)
- for each design, we run NRA for each model and store the value of the predicted pca yield wrt glucose uptake
- in this way we get a matrix of solutions. for example if we have 30 models and 300 unique designs across all models,
we will have a matrix of 30x300 solutions.
"""
import numpy as np
import pandas as pd
import glob

from pytfa.io.json import load_json_model
from pytfa.optim.constraints import *

from skimpy.io.yaml import load_yaml_model
from skimpy.io.regulation import load_enzyme_regulation
from skimpy.analysis.oracle.load_pytfa_solution import load_fluxes, load_concentrations, load_equilibrium_constants
from skimpy.analysis.mca.utils import get_dep_indep_vars_from_basis
from skimpy.utils.general import get_stoichiometry
from skimpy.core.parameters import load_parameter_population
from skimpy.utils.tabdict import TabDict
from skimpy.utils.namespace import NET, QSSA

from NRAplus.core.model import NRAmodel

from sympy import Matrix
from scipy.sparse import csc_matrix as sparse_matrix


# Design parameters
N_ENZYMES = 3
MAX_FOLD_ENZ_ACTIVITY = 2
MAX_FOLD_CONC_CHANGE = 5

# Cellular parameters
CONCENTRATION_SCALING = 1e9 # 1 mol to 1 mmol
TIME_SCALING = 1 # 1hr
DENSITY = 1200 # g/L
GDW_GWW_RATIO = 0.3 # Assumes 70% Water
flux_scaling_factor = 1e-3 * (GDW_GWW_RATIO * DENSITY) * CONCENTRATION_SCALING / TIME_SCALING

# Paths
path_to_tmodel = './../../models/constraint_based_model_scaffold.json'
path_to_kmodel = './../../models/kinetic_model_scaffold.yaml'
path_to_tfa_sample = './../../data/steady_states.csv'
path_to_designs = '../output/designs/'
path_to_conservation_relations = '../../data/ST10284_cons_relations_annotated.csv'
path_to_params = './../../kinetic-modelling/parameters/parameters_top_9_models.hdf5'

# Load models, samples and parameters
tmodel = load_json_model(path_to_tmodel)
kmodel = load_yaml_model(path_to_kmodel)

# Get the correct parameter sets and tfa sample
tfa_sample = pd.read_csv(path_to_tfa_sample, index_col=0).iloc[3191] # steady state closes to mean
parameter_population = load_parameter_population(path_to_params)
kinetic_models = parameter_population._index.keys()

# Load all NRA designs for all models
list_designs = []
for this_model in kinetic_models:
    folder = path_to_designs + '{}'.format(this_model)
    df = pd.read_csv(folder + '/designs.csv', index_col = 0) # IMP - either designs or designs_only.csv
    list_designs.append(df)
df_designs = pd.concat(list_designs)

# Drop solutions and convert all designs to binary (1 if the enzyme is targeted in a design and 0 otherwise)
df_designs = df_designs.notnull().astype("int")
df_designs = df_designs.drop('solution', axis=1)

# Get only unique designs and convert into a list of designs
df_unique_designs = df_designs.drop_duplicates()
list_of_designs = []
for _, this_design in df_unique_designs.iterrows():
    this_design = this_design[this_design == 1]
    list_of_designs.append(list(this_design.index))

print("Number of unique designs is : {}".format(len(list_of_designs)))

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
# Prepare the kinetic models for MCA
kmodel.prepare(mca=False)

conservations = pd.read_csv(path_to_conservation_relations, index_col=0)
# Transform the conservation realtsion to relfect the dependent and independent weights
# (use set of orthogonal dependents).
L0, pivot = Matrix(conservations[kmodel.reactants].values).rref()
kmodel.conservation_relation = sparse_matrix(L0,dtype=np.float)
dep_ix, indep_ix = get_dep_indep_vars_from_basis(kmodel.conservation_relation)
kmodel.independent_variables_ix = indep_ix
kmodel.dependent_variables_ix = dep_ix
kmodel.reduced_stoichiometry = get_stoichiometry(kmodel,kmodel.reactants)[indep_ix,:]
kmodel.compile_jacobian(ncpu=4)

kmodel.compile_mca(mca_type=NET, sim_type=QSSA, parameter_list=parameter_list, ncpu=10)

df_fold_changes = [] # this list will become a dataframe for the fold changes for each combination of model and design
for this_model in kinetic_models:
    print(this_model)
    tfa_ix, kin_ix = this_model.split(',')

    # DGo of this reaction needs to be relaxed during design generation
    tmodel.variables['DGo_E4Ptm'].lb = -2

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

    # Get control coefficients - remember multiply fluxes by 3 for consistency with 3x vmax
    flux_control_coeff = kmodel.flux_control_fun(fluxes * 3, concentrations, [parameter_set])
    concentration_control_coeff = kmodel.concentration_control_fun(fluxes * 3, concentrations, [parameter_set])
    ccc = concentration_control_coeff.mean('sample')
    fcc = flux_control_coeff.mean('sample')

    def reset_nmodel(tmodel, this_tfa_sample, fcc, ccc):
        """ function for repeated calls to create nmodel
        """

        # Remove bounds on LCs concentrations in the tmodel first
        for this_LC in tmodel.log_concentration:
            this_LC.variable.ub = 100
            this_LC.variable.lb = -100

        # Build NRAmodel based on the mean control coeffcients
        nmodel = NRAmodel(tmodel, this_tfa_sample, fcc, ccc, type=NET)

        nmodel.objective = nmodel.variables.LFR_PCAtex - nmodel.variables.LFR_GLCt1  # anthranilate yield change
        nmodel.variables.LFR_LMPD_s_0450_c_1_256.lb = -0.2

        # General Design constraints
        nmodel.max_enzyme_modifications = N_ENZYMES
        nmodel.max_fold_enzyme_change = np.log(MAX_FOLD_ENZ_ACTIVITY)
        nmodel.max_fold_concentration_change = np.log(MAX_FOLD_CONC_CHANGE)
        for var in nmodel.enzyme_down_regulation_variable:  # Allow for knockouts
            var.ub = 100
        nmodel.variables.EU_DDPAm.ub = 0
        nmodel.variables.ED_DDPAm.ub = 0

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
    Test each unique design for this kinetic model
    """

    this_nmodel = nmodel.copy()
    for ix, this_design in enumerate(list_of_designs):

        # Enforce current design (enzymes on/off)
        for this_enzyme in this_design:
            this_nmodel.variables[this_enzyme].lb = 1e-3 # Enforcing this design
        try:
            sol = this_nmodel.optimize()
            obj = sol.objective_value
        except:
            print("ix {} did not work".format(ix))
            obj = -100

        # Reset enzyme regulations
        for this_enzyme in this_design:
            this_nmodel.variables[this_enzyme].lb = 0

        # Store all the design fold changes for later use (for bioreactor simulations)
        for this_enzyme in this_design:
            dict_ = {'model': this_model,
                     'design_ix': ix,
                     'enzyme': this_enzyme,
                     'fold-change': sol.raw.loc[this_enzyme],
                     'nra_sol': obj,
                     }
            df_fold_changes.append(dict_)

# Create a df of all nra solutions - index = all the unique designs | columns = each kinetic model
df_fold_changes = pd.DataFrame(df_fold_changes)
df_fold_changes.to_csv(path_to_designs + '/nra_sol_all_designs_all_models.csv')