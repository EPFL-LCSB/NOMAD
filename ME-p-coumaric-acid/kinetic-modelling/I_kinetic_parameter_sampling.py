'''
In this script, we generate kinetic parameter sets around the TFA sample thta is closest to the mean of all TFA samples (3191)
KM bounds are integrated from BRENDA with +- 1 order of magnitude 
i.e. lowest KM over all the yeast KMs *0.1 and highest KM over all yeast KMs *10
'''

from pytfa.io.json import load_json_model

from skimpy.io.yaml import load_yaml_model
from skimpy.analysis.oracle.load_pytfa_solution import load_fluxes, \
    load_concentrations, load_equilibrium_constants
from skimpy.core.parameters import ParameterValuePopulation
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


NCPU = 2
N_SAMPLES = 100 # Exploitation of the TFA sample closest to the mean - 3191

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
path_to_samples = './../data/steady_states.csv'
path_for_output = './parameters/parameters_sample_3191_.hdf5' # Adding an _ at the end so that yo don't overwrite the models we have generated

# Load pytfa model
tmodel = load_json_model(path_to_tmodel)

# Load the kinetic model
kmodel = load_yaml_model(path_to_kmodel)

# Prep and compile
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

kmodel.compile_jacobian(ncpu=NCPU)

# Load tfa samples
samples = pd.read_csv(path_to_samples, header=0, index_col=0).iloc[3191:3192]

# Initialize parameter sampler
sampler_params = SimpleParameterSampler.Parameters(n_samples=N_SAMPLES)
sampler = SimpleParameterSampler(sampler_params)

# Set bounds on Kms
kms = [p for p in kmodel.parameters.keys() if p.startswith('km_')]
for this_km in kms:
    kmodel.parameters[this_km].bounds = (KM_LB * CONCENTRATION_SCALING, KM_UB * CONCENTRATION_SCALING)

S = get_stoichiometry(kmodel, kmodel.reactants).todense()

lambda_max_all = []
lambda_min_all = []
good_samples = []
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
                                                  max_trials=15000,
                                                  bounds_sample=(0.001, 0.999), # bias the bounds away from low sat
                                                  )
    
    lambda_max_all.append(pd.DataFrame(lamda_max))
    lambda_min_all.append(pd.DataFrame(lamda_min))

    # Get good models only i.e. those that are 3x faster than the doubling time
    good_indices = [i for i in range(len(lamda_max)) if lamda_max[i] < -0.4]
    if len(good_indices) > 0:
        params_population = ParameterValuePopulation([params[i] for i in good_indices], kmodel)
        params_population.save(path_for_output.format(i))

# Process df and save dataframe
lambda_max_all = pd.concat(lambda_max_all, axis=1)
lambda_min_all = pd.concat(lambda_min_all, axis=1)

# lambda_max_all.to_csv('./eigenvalues/maximal_eigenvalues_3191.csv'.format(LOWER_IX,UPPER_IX))
# lambda_min_all.to_csv('./eigenvalues/minimal_eigenvalues_3191.csv'.format(LOWER_IX,UPPER_IX))







