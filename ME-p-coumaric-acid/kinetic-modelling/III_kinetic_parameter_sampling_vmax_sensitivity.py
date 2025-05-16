'''
Vmax sensitivity analysis
1. 3x all fluxes.
2. Fix all concs + thermo vars to the 3191 steady state (old one) 
3. Fix exofluxomics
4. Allow 20% change in intracellular fluxes, concs and thermo vars
5. Sample fluxes, concs, thermo vars 
6. Build models around these using same Kms as for one of the the initial 9 models (3191,8)
7. Check distribution of Vmaxs

Here we do step 6 for the distribution of Vmaxs
'''

from pytfa.io.json import load_json_model

from skimpy.io.yaml import load_yaml_model
from skimpy.analysis.oracle.load_pytfa_solution import load_fluxes, \
    load_concentrations, load_equilibrium_constants
from skimpy.core.parameters import ParameterValuePopulation
from skimpy.core.parameters import load_parameter_population
from skimpy.sampling.simple_parameter_sampler import SimpleParameterSampler
from skimpy.sampling.ga_parameter_sampler import GaParameterSampler
from skimpy.utils.general import get_stoichiometry
from skimpy.io.regulation import load_enzyme_regulation
from skimpy.analysis.mca.utils import get_dep_indep_vars_from_basis

from sympy import Matrix
from scipy.sparse import csc_matrix as sparse_matrix

import pandas as pd
import numpy as np
from sys import argv
import os


# We will generate models for each of the 100 steady-states
LOWER_IX = 0
UPPER_IX = 100

NCPU = 4
N_SAMPLES = 20 # 20 kinetic parameters sets per steady-state

CONCENTRATION_SCALING = 1e9 # 1 mol to 1 mmol
TIME_SCALING = 1 # 1hr to 1min

# Parameters of the Yeast cell
DENSITY = 1200 # g/L
GDW_GWW_RATIO = 0.3 # Assumes 70% Water
KM_LB = 3.5e-10 # min. for S cerevisiae in BRENDA
KM_UB = 5 # max for S.cerevisiae in BRENDA
flux_scaling_factor = 1e-3 * (GDW_GWW_RATIO * DENSITY) * CONCENTRATION_SCALING / TIME_SCALING # mmol/gDW/hr

# Paths
path_to_kmodel = './../models/kinetic_model_scaffold.yaml'
path_to_tmodel = './../models/constraint_based_model_scaffold.json'
path_to_samples = './../data/steady_states_vmax_sensitivity.csv'
path_to_old_param_set = './parameters/parameters_top_9_models.hdf5'
path_for_output = './parameters/vmax-sensitivity-/' # adding the - at the end to ensure no overwriting

if not os.path.exists(path_for_output):
    os.makedirs(path_for_output)

# Load pytfa model
tmodel = load_json_model(path_to_tmodel)

# Load the kinetic model
kmodel = load_yaml_model(path_to_kmodel)
parameter_population = load_parameter_population(path_to_old_param_set)
old_param_set = parameter_population['3191,8'] # 3191,8 to keep the KMs from this model

# Prep and compile
# Manual matlab curation
kmodel.prepare(mca=False)

conservations = pd.read_csv('./../data/ST10284_cons_relations_annotated.csv',index_col=0)

# Transform the conservation realtsion to relfect the dependent and independent weights
# (use set of orthogonal dependents).
L0, pivot = Matrix(conservations[kmodel.reactants].values).rref()
kmodel.conservation_relation = sparse_matrix(L0,dtype=np.float)
dep_ix, indep_ix = get_dep_indep_vars_from_basis(kmodel.conservation_relation)
kmodel.independent_variables_ix = indep_ix
kmodel.dependent_variables_ix = dep_ix
kmodel.reduced_stoichiometry = get_stoichiometry(kmodel,kmodel.reactants)[indep_ix,:]

kmodel.compile_jacobian(ncpu=10)

# Load tfa samples
samples = pd.read_csv(path_to_samples, header=0, index_col=0)

# Initialize parameter sampler
sampler_params = SimpleParameterSampler.Parameters(n_samples=N_SAMPLES)
sampler = SimpleParameterSampler(sampler_params)

# Set bounds on Kms from one of the old models
kms = [p for p in kmodel.parameters.keys() if p.startswith('km_')]
for this_km in kms:
    kmodel.parameters[this_km].bounds = (old_param_set[this_km] * 0.99, old_param_set[this_km] * 1.01)

S = get_stoichiometry(kmodel, kmodel.reactants).todense()

for i, sample in samples.iterrows():

    # Load fluxes and concentrations
    fluxes = load_fluxes(sample, tmodel, kmodel,
                         density=DENSITY,
                         ratio_gdw_gww=GDW_GWW_RATIO,
                         concentration_scaling=CONCENTRATION_SCALING,
                         time_scaling=TIME_SCALING)

    concentrations = load_concentrations(sample, tmodel, kmodel,
                                         concentration_scaling=CONCENTRATION_SCALING)

    # Fetch equilibrium constants
    load_equilibrium_constants(sample, tmodel, kmodel,
                               concentration_scaling=CONCENTRATION_SCALING,
                               in_place=True)

    params, lamda_max, lamda_min = sampler.sample(kmodel, fluxes, concentrations,
                                                  only_stable=True,
                                                  min_max_eigenvalues=True,
                                                  max_trials=1000,
                                                  bounds_sample=(0.001, 0.999), # bias the bounds away from low sat
                                                  )
    
    try:
        params_population = ParameterValuePopulation(params, kmodel)
        params_population.save(path_for_output + 'parameters_sample_{}.hdf5'.format(i))

    except:
        print('Could not create parameters for this sample!')
        continue






