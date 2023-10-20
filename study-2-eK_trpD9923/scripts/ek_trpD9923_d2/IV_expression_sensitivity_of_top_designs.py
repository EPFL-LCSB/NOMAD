"""
In this script we conduct and enzyme expression sensitivity analysis of the top 5 designs in all the 10 models in a nonlinear batch setting.
We test a +- 50% perturabation to all 3 enzymes simultaneously around the mean NRA suggested values
This test sufficed as all the designs passed it.

"""

from pytfa.io.json import load_json_model
from skimpy.io.yaml import load_yaml_model

from skimpy.simulations.reactor import make_batch_reactor
from skimpy.core.solution import ODESolutionPopulation

from skimpy.analysis.oracle.load_pytfa_solution import load_fluxes, load_concentrations, load_equilibrium_constants
from skimpy.core.parameters import load_parameter_population
from skimpy.utils.tabdict import TabDict
from skimpy.utils.namespace import NET, QSSA
from skimpy.io.regulation import load_enzyme_regulation
from skimpy.analysis.ode.utils import make_flux_fun

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
from sys import argv
import os


'''
Sensitivty analysis parameters
- 0 = EU_ANS
- 1 = EU_DDPA / EU_DHQS
- 2 = third enzymes -> the other one
- 3 = all three enzymes
'''
ENZ_TO_PERTURB = 3
PERCENT_PERTURBATION = 50
TOTAL_TIME = 60
N_STEPS = 1000
MAX_FOLD_CONC_CHANGE = 2.5

'''
Paths
'''
path_to_tmodel = '../../models/tfa_model.json'
path_to_kmodel = '../../models/kinetic_model_scaffold.yaml'
path_to_tfa_samples = '../../data/steady_state_samples.csv'
path_to_params = '../../data/enhanced_kinetic_models.hdf5'
path_to_regulation_data = '../../data/allosteric_regulations.csv'
folder_to_designs = './../../output/data/eK_trpD9923_d2'
path_to_designs = folder_to_designs + '/all_unique_designs_w_enz_fold_changes.csv'

if ENZ_TO_PERTURB < 3:
    folder_for_output = './../../output/enzyme-sensitivity/eK_trpD9923_d2/{}-percent/enzyme-{}'.format(PERCENT_PERTURBATION, ENZ_TO_PERTURB)
else:
    folder_for_output = './../../output/enzyme-sensitivity/eK_trpD9923_d2/{}-percent/all-enzymes'.format(PERCENT_PERTURBATION)

if not os.path.exists(folder_for_output):
    os.makedirs(folder_for_output)

# Cellular parameters
CONCENTRATION_SCALING = 1e9 # 1 mol to 1 mmol
TIME_SCALING = 1 # 1hr
DENSITY = 1105 # g/L
GDW_GWW_RATIO = 0.3 # Assumes 70% Water
flux_scaling_factor = 1e-3 * (GDW_GWW_RATIO * DENSITY) * CONCENTRATION_SCALING / TIME_SCALING

'''
Model and simulation initialization and preparation
'''
# Load tfa models, samples, and parameters
tmodel = load_json_model(path_to_tmodel)
kmodel_draft = load_yaml_model(path_to_kmodel)
samples = pd.read_csv(path_to_tfa_samples, header=0, index_col=0)

# Load the dataframe with regulations and choose those rxns and mets that are in the kinetic model
df = pd.read_csv(path_to_regulation_data)
df_reg = df[df['reaction_id'].isin(list(kmodel_draft.reactions.keys()))]
df_reg = df_reg[df_reg['regulator'].isin(list(kmodel_draft.reactants.keys()))]

kmodel = load_enzyme_regulation(kmodel_draft, df_reg)

reactor = make_batch_reactor('single_species.yaml', df_regulation= df_reg)
reactor.compile_ode(add_dilution=False)
reactor_volume = reactor.models.strain_1.parameters.strain_1_volume_e.value

# Get the correct parameter sets
parameter_population = load_parameter_population(path_to_params)
kinetic_models = list(parameter_population._index.keys())

list_all = []
list_of_flux_dfs = []
dict_number_good = {}

# Get the top designs and fold changes
# Identify top designs - ones with the best mean
df_designs = pd.read_csv(path_to_designs, header=0, index_col=0)
top_design_indices = df_designs.groupby('design_ix').mean().sort_values(by='nra_sol', ascending=False).index[:5]
top_designs = df_designs[df_designs['design_ix'].isin(top_design_indices)]

'''
Sensitivity analysis
- Go through each design. For each design go through the models, generate a vector of vmax perturbations for the enzymes
- Simulate the perturbed system
- Store the solutions
'''
for design_ix, this_design in top_designs.groupby('design_ix'):
    
    # Get mean NRA fold changes for this design across the 10 models
    mean_fold_changes = this_design.groupby('enzyme').mean()
    this_design_sols = []
    this_design_indices = []
    df_fold_changes = []

    # Loop through all the models and for each model simulate 10 perturbed systems
    for ix in kinetic_models:
        # Get reference concentrations
        tfa_id,_ = ix.split(',')
        tfa_id = int(tfa_id)
        sample = samples.loc[tfa_id]
        reference_concentrations = load_concentrations(sample, tmodel, kmodel, concentration_scaling=reactor.concentration_scaling)

        # Load fluxes and concentrations
        fluxes = load_fluxes(sample, tmodel, kmodel,
                             density=DENSITY,
                             ratio_gdw_gww=GDW_GWW_RATIO,
                             concentration_scaling=CONCENTRATION_SCALING,
                             time_scaling=TIME_SCALING)

        # Fetch equilibrium constants
        load_equilibrium_constants(sample, tmodel, kmodel,
                                   concentration_scaling=CONCENTRATION_SCALING,
                                   in_place=True)


        def reset_reactor(reactor, reference_concentrations, parameter_population, sample,):

            reactor.parametrize(parameter_population[ix], 'strain_1')
            reactor.initialize(reference_concentrations, 'strain_1')

            # Initial biomass : from the paper it seemed to be 0.037 g/L - wt of a cell is 0.28e-12 g
            reactor.initial_conditions['biomass_strain_1'] = 0.037 *0.05 / 0.28e-12

            for met_ in reactor.medium.keys():
                LC_id = 'LC_' + met_
                LC_id = LC_id.replace('_L','-L')
                LC_id = LC_id.replace('_D', '-D')
                reactor.initial_conditions[met_] = np.exp(sample.loc[LC_id]) * 1e9

            # Volume parameters
            reactor.models.strain_1.parameters.strain_1_volume_e.value = reactor_volume
            reactor.models.strain_1.parameters.strain_1_cell_volume_e.value = 1.0 # 1.0 #(mum**3) look up cell volume bionumbers
            reactor.models.strain_1.parameters.strain_1_cell_volume_c.value = 1.0 # 1.0 #(mum**3)
            reactor.models.strain_1.parameters.strain_1_cell_volume_p.value = 1.0  # 1.0 #(mum**3)
            reactor.models.strain_1.parameters.strain_1_volume_c.value = 0.9*1.0 #(mum**3)
            reactor.models.strain_1.parameters.strain_1_volume_p.value = 0.1*1.0 # (mum**3)

            return reactor

        def print_yields(sol):

            time = np.linspace(0, TOTAL_TIME, N_STEPS)
            biomass = sol.concentrations['biomass_strain_1'] * 0.28e-12
            glc = sol.concentrations['glc_D_e'] * 180.156 * 1e-9
            ix_end_ferm = np.where(biomass >= 0.99 * np.max(biomass))[0][0]
            glc_consumed = glc[0] - glc[ix_end_ferm]
            tot_glc_consumed = glc[0] - glc.iloc[-1]
            final_biomass = sol.concentrations.loc[len(sol.concentrations) - 1, 'biomass_strain_1'] * 0.28e-12 / 0.05
            final_anth_conc = sol.concentrations.loc[len(sol.concentrations) - 1, 'anth_e'] * 136.13 * 1e-9

            print('{}: biomass = {}, ferm time = {}, tot glc consumed = {}, anth_titer = {}'.format(ix,
                                                                                                    final_biomass,
                                                                                                    time[ix_end_ferm],
                                                                                                    tot_glc_consumed,
                                                                                                    final_anth_conc))

        """
        Get wildtype solution first
        """
        reactor = reset_reactor(reactor, reference_concentrations, parameter_population, sample)
        sol_wt = reactor.solve_ode(np.linspace(0, TOTAL_TIME, N_STEPS),
                                  solver_type='cvode',
                                  rtol=1e-9,
                                  atol=1e-9,
                                  max_steps=1e9,
                                  )

        """
        Simulate perturbations to the system for this design and this model
        """
        import time

        for j in range(10):

            reactor = reset_reactor(reactor, reference_concentrations, parameter_population, sample)

            # Change the DDPA inhibition/tktA upregulation to correspond to eK_trpD9923_d2
            reactor.parameters['strain_1_ki_inhibitor1_DDPA'].value = 1e42
            reactor.parameters['strain_1_k_inhibition_IM_phe_L_c_DDPA'].value = 1e42
            reactor.parameters['strain_1_vmax_forward_TKT1'].value *= 10
            reactor.parameters['strain_1_vmax_forward_TKT2'].value *= 10
            
            # Create this temporary list to ensure that ANS and DDPA are enzymes 1 and 2
            if 'EU_DDPA' in mean_fold_changes.index:
                temp_enzyme_list = ['EU_ANS', 'EU_DDPA']
                reactor.parameters['strain_1_ki_inhibitor1_DDPA'].value = 1e42
                reactor.parameters['strain_1_k_inhibition_IM_phe_L_c_DDPA'].value = 1e42
            elif 'EU_DHQS' in mean_fold_changes.index:
                temp_enzyme_list = ['EU_ANS', 'EU_DHQS']
            temp_enzyme_list.extend([i for i in mean_fold_changes.index if i not in temp_enzyme_list])
            for this_enz in mean_fold_changes.index:
                # If ENZ_TO_CHECK <3  then you perturb that particular enzyme, otherwise perturb all 3
                if ENZ_TO_PERTURB == 3:
                    perturbation = 1 + (np.random.uniform(-1, 1) * PERCENT_PERTURBATION / 100)
                elif this_enz == temp_enzyme_list[ENZ_TO_PERTURB]:
                    perturbation = 1 + (np.random.uniform(-1, 1) * PERCENT_PERTURBATION / 100)
                else:
                    perturbation = 1

                fold_change =  perturbation * np.exp(mean_fold_changes.loc[this_enz, 'fold-change'])
                if this_enz.startswith('EU_'):
                    reactor.parameters['strain_1_vmax_forward_' + this_enz[3:]].value *= fold_change
                if this_enz.startswith('ED_'):
                    reactor.parameters['strain_1_vmax_forward_' + this_enz[3:]].value *= 1/fold_change

                # Store the details of the enzyme that was perturbed and the fold perturbation applied
                dict_info = {'model': ix,
                             'perturbation': j,
                             'design': design_ix,
                             this_enz + 'mean': np.exp(mean_fold_changes.loc[this_enz, 'fold-change']),
                             this_enz + 'perturbation': perturbation,
                             }
                df_fold_changes.append(dict_info)

            # Solve the ODEs for the current perturbed system
            start = time.time()

            if hasattr(reactor, 'solver'):
                delattr(reactor, 'solver')

            sol_1 = reactor.solve_ode(np.linspace(0, TOTAL_TIME, N_STEPS),
                                    solver_type='cvode',
                                    rtol=1e-9,
                                    atol=1e-9,
                                    max_steps=1e9,
                                    )

            end = time.time()
            print_yields(sol_1)
            this_design_sols.append(sol_1)
            this_design_indices.append(ix + ',' + str(j))

    solpop = ODESolutionPopulation(this_design_sols, this_design_indices)
    solpop.data.to_csv(folder_for_output + '/ode_sols_design_{}.csv'.format(design_ix))
    del solpop

    df_fold_changes = pd.DataFrame(df_fold_changes)
    df_fold_changes.to_csv(folder_for_output + '/fold_changes_design_{}.csv'.format(design_ix))

