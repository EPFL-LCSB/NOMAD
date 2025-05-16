"""
This script extracts the relevant data for analysis
The info from this is used for Figure 2A
1. pca titer and yield wrt glucose at 5 hours

"""
import pandas as pd
import numpy as np


CONCENTRATION_SCALING = 1e9 # 1 mol to 1 mmol
DENSITY = 1200 # g/L
GDW_GWW_RATIO = 0.3 # Assumes 70% Water
CELL_VOLUME = 42 # mum3
REACTOR_VOLUME = 0.0025 #L
GDW_PER_CELL = GDW_GWW_RATIO * DENSITY * CELL_VOLUME * 1e-15 # 1e-15 L per mum3
TOTAL_TIME = 48


# Paths
folder_to_sols = './../output/nonlinear-verification/'
path_to_designs = './../output/designs/nra_sol_all_designs_all_models.csv'

kinetic_models = ['3191,8', '3191,42', '3191,60', '3191,61', '3191,95', '3191,99',
                   '3191,104', '3191,139', '3191,142']


pca_scaling = 1e-9 * 164 # mwt of pca is 164 and concentration scaling for the models is 1e9
glc_scaling = 1e-9 * 180 # m wt of glucse is 180

# Load data from NRA
all_designs = pd.read_csv(path_to_designs, header=0, index_col=0)
top_design_indices = all_designs.groupby('design_ix').mean().sort_values(by='nra_sol', ascending=False).index[:10]
top_designs_all_models = all_designs[all_designs['design_ix'].isin(top_design_indices)]
top_designs_nra_sols = top_designs_all_models.groupby('design_ix').nra_sol.mean()
top_designs = top_designs_all_models[top_designs_all_models.model == '3191,8'][['design_ix', 'enzyme', 'nra_sol']]

# Pivot it to convert from long to wide format
# Assigning a sequence number for each enzyme within its group by design_ix
top_designs['enzyme_no'] = top_designs.groupby('design_ix').cumcount() + 1

# Pivoting the DataFrame
top_designs = top_designs.pivot(index='design_ix', columns='enzyme_no', values='enzyme')

def get_perf_params(t=24):
    '''
    function to get the performance parameters at time t
    '''
    df_ = []
    top_designs_ = top_designs.copy()

    # Loop through all models
    for m in kinetic_models:
        
        folder = folder_to_sols + '/{}'.format(m)
        
        # Load solutions
        wt = pd.read_csv(folder + '/sol_wt.csv')
        designs = pd.read_csv(folder + '/sol_designs.csv', index_col=0) # recombinant strains
        
        # Get spec prod, and yield on glucose at time t, along with pca titer for the reference strain!
        glc_wt_t0 = wt.iloc[0]['glc__D_e']
        ix = np.where(designs[designs.solution_id == 0].time >= t)[0][0] # Index where time = t
        wt_sol = wt.iloc[ix]
        # Remember biomass is in number of cells, need to convert to gDW for specific productivity
        spec_prod_wt = (wt_sol.loc['pca_e']/ CONCENTRATION_SCALING) / (wt_sol.loc['biomass_strain_1'] * GDW_PER_CELL) / t
        glc_yield_wt = (wt_sol.loc['pca_e'] * pca_scaling) / ((glc_wt_t0 - wt_sol.loc['glc__D_e']) * glc_scaling)
        des_sols = []
        
        # Do the same for the recombinant straints after checking if they converged
        for ix_design, design in designs.groupby('solution_id'):
            converged = 0
            if design.time.iloc[-1] == TOTAL_TIME:
                converged = 1
            else:
                print('{},{} has not converged!!'.format(m, ix_design))
                continue

            des_sol = design.iloc[ix]
            spec_prod_des = (des_sol.loc['pca_e']/ CONCENTRATION_SCALING) / (des_sol.loc['biomass_strain_1'] * GDW_PER_CELL) / t
            glc_yield_des = (des_sol.loc['pca_e'] * pca_scaling) / ((glc_wt_t0 - des_sol.loc['glc__D_e']) * glc_scaling)
            des_sol = des_sol[['solution_id', 'biomass_strain_1', 'pca_e']].copy()
            des_sol['spec_prod'] = spec_prod_des
            des_sol['glc_yield'] = glc_yield_des
            des_sol['converged'] = converged
            des_sols.append(des_sol)            
        
        # Now get the increase in pca yield on glucose, specific productivity and pca titer wrt wildtype at that time.
        des_sols = pd.concat(des_sols, axis=1).T
        des_sols['pca_inc'] = des_sols['pca_e'] / wt_sol.loc['pca_e']
        des_sols['specific_prod_inc'] = des_sols['spec_prod'] / spec_prod_wt
        des_sols['glc_yield_inc'] = des_sols['glc_yield'] / glc_yield_wt
        des_sols['model'] = m
        
        df_.append(des_sols)
        
    df_designs = pd.concat(df_)
    df_designs.columns = ['design', 'biomass', 'pca', 'specific_prod', 'glc_yield', 'converged', 'pca_inc', 'specific_prod_inc', 'glc_yield_inc', 'model']
    df_designs.design = df_designs.design.astype('int')

    # How many designs give an inc in pca specific productivity?
    inc_sprod = df_designs[df_designs.specific_prod_inc.gt(1)]
    inc_glc_yield = df_designs[df_designs.glc_yield_inc.gt(1)]
    inc_pca = df_designs[df_designs.pca_inc.gt(1)]

    # Collate design info (enzyme membership) and the performance in nonlin simulations
    top_designs_['perf_pca'] = inc_pca.design.value_counts()
    top_designs_['perf_glc_yield'] = inc_glc_yield.design.value_counts()
    top_designs_['nra_sol'] = top_designs_nra_sols
    top_designs_['avg_pca_inc'] = df_designs.groupby('design').pca_inc.mean()
    top_designs_['avg_glc_yield_inc'] = df_designs.groupby('design').glc_yield_inc.mean()
    top_designs_['strain'] = pd.Series(comp_exp_links)

    return df_designs, top_designs_

# Get info at 5 hours and 20 hours for example
df_designs_5h, top_designs_5h = get_perf_params(t=5)
df_designs_20h, top_designs_20h = get_perf_params(t=20)

df_designs_5h.to_csv('./../output/nonlinear-verification/performance_analysis.csv')
