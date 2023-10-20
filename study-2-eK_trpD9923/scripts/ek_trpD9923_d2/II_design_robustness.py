"""
In this script, we conduct design robustness analysis **for the d2 design**
- we take in the kinetic models, and all NRA designs across all models
- we then prune for unique designs (unique by membership, not value)
- for each design, we run the NRA tool for each model and store the value
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
from skimpy.core.parameters import load_parameter_population
from skimpy.utils.tabdict import TabDict
from skimpy.utils.namespace import NET, QSSA

from NRAplus.core.model import NRAmodel


# Design parameters
N_ENZYMES = 3
MAX_FOLD_ENZ_ACTIVITY = 10
MAX_FOLD_CONC_CHANGE = 2.5

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
folder_to_designs = './../../output/data/eK_trpD9923_d2'

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
list_designs = []
for this_model in kinetic_models:
    output_folder = folder_to_designs + '/{}'.format(this_model)
    df = pd.read_csv(output_folder + '/designs.csv', index_col = 0) # IMP - either designs or designs_only.csv
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
kmodel.prepare()
kmodel.compile_mca(mca_type=NET, sim_type=QSSA, parameter_list=parameter_list, ncpu=4)

dict_all = {} # Dictionary to store all nra solutions obtained using each design, for each kmodel
df_fold_changes = []
for this_model in kinetic_models:
    print(this_model)
    tfa_ix, _ = this_model.split(',')

    this_tfa_sample = pd.read_csv(path_to_tfa_samples, header=0, index_col=0).iloc[int(tfa_ix)]
    parameter_set = parameter_population['{}'.format(this_model)]
    kmodel.parameters = parameter_set

    # Load fluxes and concentrations of the double mutant as pandas series - this is important
    fluxes = pd.read_csv(folder_to_designs + '/{}/steady_state_fluxes.csv'.format(this_model), index_col=0).squeeze()
    concentrations = pd.read_csv(folder_to_designs + '/{}/steady_state_concentrations.csv'.format(this_model), index_col=0).squeeze()
    new_parameter_set = {kmodel.parameters[k].symbol:v for k,v in parameter_set.items()}
    new_parameter_set[kmodel.parameters['ki_inhibitor1_DDPA'].symbol] = 1e42
    new_parameter_set[kmodel.parameters['k_inhibition_IM_phe_L_c_DDPA'].symbol] = 1e42
    new_parameter_set[kmodel.parameters['vmax_forward_TKT1'].symbol] *= 10
    new_parameter_set[kmodel.parameters['vmax_forward_TKT2'].symbol] *= 10

    # Fetch equilibrium constants
    load_equilibrium_constants(this_tfa_sample, tmodel, kmodel,
                               concentration_scaling=CONCENTRATION_SCALING,
                               in_place=True)

    # Get control coefficients
    flux_control_coeff = kmodel.flux_control_fun(fluxes, concentrations, [new_parameter_set])
    concentration_control_coeff = kmodel.concentration_control_fun(fluxes, concentrations, [new_parameter_set])
    ccc = concentration_control_coeff.mean('sample')
    fcc = flux_control_coeff.mean('sample')

    def reset_nmodel(tmodel, this_tfa_sample, fcc, ccc):
        """ function for repeated calls to create nmodel
        """

        # Remove bounds on anthranilate concentrations in the tmodel first
        for this_LC in tmodel.log_concentration:
            this_LC.variable.ub = 100
            this_LC.variable.lb = -100

        # Build NRAmodel based on the mean control coeffcients
        nmodel = NRAmodel(tmodel, this_tfa_sample, fcc, ccc, type=NET)

        # Set Objective maximise the relative change of anthranilate secretion
        nmodel.objective = nmodel.variables.LFR_ANS - nmodel.variables.LFR_GLCtex

        # Ensure that TKT is expressed together
        expr = nmodel.variables.EU_TKT1 - nmodel.variables.EU_TKT2
        nmodel.add_constraint(kind=ModelConstraint, id_='TKT_ratio_up', hook=nmodel, expr=expr, lb=0, ub=0)
        expr = nmodel.variables.ED_TKT1 - nmodel.variables.ED_TKT2
        nmodel.add_constraint(kind=ModelConstraint, id_ = 'TKT_ratio_down', hook=nmodel, expr= expr, lb = 0, ub = 0)

        # Design constraints
        nmodel.max_enzyme_modifications = N_ENZYMES
        nmodel.max_fold_enzyme_change = np.log(MAX_FOLD_ENZ_ACTIVITY)
        nmodel.max_fold_concentration_change = np.log(MAX_FOLD_CONC_CHANGE)
        for var in nmodel.enzyme_down_regulation_variable:  # Allow for knockouts
            var.ub = 100
        nmodel.variables.EU_LMPD_biomass_c_1_420.ub = 1e-3  # Don't allow enzymatic changes to the lumped reaction
        nmodel.variables.ED_LMPD_biomass_c_1_420.ub = 1e-3

        # Numerical configurations
        nmodel.solver.configuration.tolerances.feasibility = 1e-9
        nmodel.solver.configuration.tolerances.optimality = 1e-9
        nmodel.solver.configuration.tolerances.integrality = 1e-9
        nmodel.solver.problem.parameters.read.scale = -1
        nmodel.solver.problem.parameters.emphasis.numerical = 1

        nmodel.solver.configuration.presolve = True #*************** VVVIMP ******************

        return nmodel

    nmodel = reset_nmodel(tmodel, this_tfa_sample, fcc, ccc)

    """
    Test each unique design for this kinetic model
    """
    dict_ = {}  # Dictionary to store each NRA solution obtained using each unique design for this kinetic model
    this_nmodel = nmodel.copy()
    for ix, this_design in enumerate(list_of_designs):

        # Enforce current design (enzymes on/off)
        for this_enzyme in this_design:
            this_nmodel.variables[this_enzyme].lb = 1e-3 # Set to tolerance level
        try:
            sol = this_nmodel.optimize()
            obj = sol.objective_value
        except:
            print("ix {} did not work".format(ix))
            obj = -100

        dict_[str(ix)] = obj
        dict_all[this_model] = dict_
        

        # Reset enzyme regulations
        for this_enzyme in this_design:
            this_nmodel.variables[this_enzyme].lb = 0

        # Store all the design fold changes for later use (for bioreactor simulations)
        for this_enzyme in this_design:
            dict__ = {'model': this_model,
                     'design_ix': ix,
                     'enzyme': this_enzyme,
                     'fold-change': sol.raw.loc[this_enzyme],
                     'nra_sol': obj,
                     }
            df_fold_changes.append(dict__)
            
# Create a df of all nra solutions - index = all the unique designs | columns = each kinetic model
df_fold_changes = pd.DataFrame(df_fold_changes)
df_all_nra_sols = pd.DataFrame(dict_all)
df_all_nra_sols.to_csv(folder_to_designs + '/all_nra_sols_across_models.csv')
df_unique_designs.to_csv(folder_to_designs + '/all_unique_designs.csv')
df_fold_changes.to_csv(folder_to_designs + '/all_unique_designs_w_enz_fold_changes.csv')

