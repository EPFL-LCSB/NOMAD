'''
This is a script to get the variation in different parameters across the 9 models.
Requested by reviewer 1.

I also take in all the BRENDA Kms, chck their coverage and see where our KMs lie in relation
'''
import numpy as np
import pandas as pd
from skimpy.io.yaml import load_yaml_model
import glob

from pytfa.io.json import load_json_model
from pytfa.optim.constraints import *

from skimpy.io.yaml import load_yaml_model

from skimpy.core.parameters import load_parameter_population
from skimpy.utils.tabdict import TabDict



# Cellular parameters
CONCENTRATION_SCALING = 1e9 # 1 mol to 1 mmol
TIME_SCALING = 1 # 1hr
DENSITY = 1200 # g/L
GDW_GWW_RATIO = 0.3 # Assumes 70% Water
CELL_VOLUME = 42 # mum3
GDW_PER_CELL = GDW_GWW_RATIO * DENSITY * CELL_VOLUME * 1e-15 # 1e-15 L per mum3
flux_scaling_factor = 1e-3 * (GDW_GWW_RATIO * DENSITY) * CONCENTRATION_SCALING / TIME_SCALING

# Paths
path_to_params = './parameters/parameters_top_9_models.hdf5'
path_to_kmodel = './../models/kinetic_model_scaffold.yaml'
path_to_brenda_kms = './../data/kms_in_our_model_manual_curation_.csv'
path_to_rxn_ec_links = './../data/reaction_ec_pairing.csv'

# Load parameters and model
kmodel = load_yaml_model(path_to_kmodel)
parameter_population = load_parameter_population(path_to_params)
kinetic_models = ['3191,8', '3191,42', '3191,60', '3191,61', '3191,95', '3191,99',
                   '3191,104', '3191,139', '3191,142']

df_vmax = []
df_km = []
for this_model in kinetic_models:
    print(this_model)
    tfa_ix, kin_ix = this_model.split(',')
    params = parameter_population[this_model]

    vmaxs = {k:v / flux_scaling_factor for k,v in params.items() if k.startswith('vmax_forward_')}
    kms = {k:v / CONCENTRATION_SCALING for k,v in params.items() if k.startswith('km_')}

    df_vmax.append(vmaxs)
    df_km.append(kms)

# Create dataframes for vmaxs and kms
df_vmax = pd.DataFrame(df_vmax)
df_vmax.index = kinetic_models

df_km = pd.DataFrame(df_km)
df_km.index = kinetic_models

# Get the description of the datasets
desc_vmax = df_vmax.agg('describe').T
desc_vmax['CV'] = desc_vmax['std'] / desc_vmax['mean']
desc_vmax.to_csv('./output/variation_in_vmaxs_.csv')
desc_km = df_km.agg('describe').T
desc_km['CV'] = desc_km['std'] / desc_km['mean']
desc_km.to_csv('./output/variation_in_kms_.csv')

'''
What is our coverage when compared to BRENDA?

'''
kms_brenda = pd.read_csv(path_to_brenda_kms, index_col=0)

# We will drop all those with empty substrate - this means I could not find a particular metabolite that is associated
kms_brenda = kms_brenda.dropna(subset=['met'])
kms_brenda.KM = kms_brenda.KM.astype(float)

# What is the coverage?
print('There are KMs for {} rxn substrate pairings across {} reactions and {} metabolites \n'.format(len(kms_brenda[['EC', 'met']].drop_duplicates()),
                                                                                                     len(kms_brenda.EC.unique()), len(kms_brenda.met.unique())))

# Initially, I was doing by EC and Substrate...but now I do EC and met.. because Substrates can be named differently
min_kms = kms_brenda.groupby(['EC', 'met']).KM.min().reset_index()
min_kms.columns = ['EC', 'met', 'KM_min']
max_kms = kms_brenda.groupby(['EC', 'met']).KM.max().reset_index()
max_kms.columns = list(['EC', 'met', 'KM_max'])
km_ratios = pd.merge(min_kms, max_kms, on=['EC', 'met'])
km_ratios['ratio'] = km_ratios['KM_max'].astype(float) / km_ratios['KM_min'].astype(float)
km_ratios.sort_values(by='ratio')
len(km_ratios[km_ratios.ratio.ge(100)])

# Where do our reaction KMs fall?
# Note: although km_ratios on;ly has 340 ish entires, kms_ has 397. This is because you can have two different reaction names with the same EC
rxn_ec_mapping = pd.read_csv(path_to_rxn_ec_links)
rxn_ec_mapping = rxn_ec_mapping[['bigg_reaction', 'eccodes']]
rxn_ec_mapping.columns = ['reaction', 'EC']
kms_ = pd.merge(rxn_ec_mapping, km_ratios, on=['EC'], how='inner')

# At this stage, we have the same EC for two different reactions, this means that some substrate pairings won't exist.
kms_['KM_model_max'] = 0
kms_['KM_model_min'] = 0

# Now for each of these reactions, go through their kmodel parameters:
# say we have EC 1.1 for ENO and EC 1.1 for FBP as well,
# Then we will have the same mets copied for both reactions but the met for ENO might not exist in FBP and vice versa..
# This is less of an issue than if have say, ATP or ADP, in both reactions - we don't know which reaction they were talking about for ATP as a substrate.
# but we will live with this.
for ix_, row in kms_.iterrows():
    # Find the corresponding substrate number
    rxn = kmodel.reactions[row.loc['reaction']]

    # We do this split '_' to account for mismatched compartments and -2 is because we know the last bit is always the compartment
    try:
        substrate_id = [k for k, v in rxn.reactants.items() if v.name.split('_')[-2] == row.loc['met']][0]
    except:
        print('Could not do {} and {}'.format(row.loc['reaction'], row.loc['met']))
        continue
    km_model_min = desc_km.loc['km_' + substrate_id + '_' + row.loc['reaction']].loc['min'] * 1e3 # Convert to mmol
    km_model_max = desc_km.loc['km_' + substrate_id + '_' + row.loc['reaction']].loc['max'] * 1e3 # Convert to mmol
    kms_.loc[ix_, 'KM_model_max'] = km_model_max
    kms_.loc[ix_, 'KM_model_min'] = km_model_min

df = kms_.copy()
# Create a new label column "rxn_Substrate"
df['label'] = df['reaction'] + '_' + df['met']

df = df[df.KM_model_max.gt(0)]

# Non overlaps
no_overlap = df[df.KM_model_min.gt(df.KM_max) | df.KM_model_max.lt(df.KM_min)]
no_overlap_gt = df[df.KM_model_min.gt(df.KM_max)]
no_overlap_lt = df[df.KM_model_max.lt(df.KM_min)]
no_overlap_gt['ratio_gt'] = no_overlap_gt['KM_model_min']  / no_overlap_gt['KM_max']
no_overlap_lt['ratio_lt'] = no_overlap_lt['KM_model_max']  / no_overlap_lt['KM_min']

# How many non overlap are within an order of magnitude?
print('Less than but within order of magnitude = {}'.format(len(no_overlap_lt[no_overlap_lt.ratio_lt.gt(0.1)])))
print('Greater than but within order of magnitude = {}'.format(len(no_overlap_gt[no_overlap_gt.ratio_gt.lt(10)])))

df.to_csv('./output/kms_brenda_vs_model_.csv')

