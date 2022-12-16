
def get_flux_ratios(nmodel,sol,reference_fluxes,flux_scaling_ratio):
    import numpy as np
    dict_fluxes_from_sol = {}
    for rxn in nmodel.reactions:
        fwd_id = 'FLF_' + rxn.id
        rev_id = 'RLF_' + rxn.id
        fwd_flux = 0
        rev_flux = 0
        try:
            fwd_flux += np.exp(sol.raw.loc[fwd_id])
        except:
            print('Forward flux does not exist for {}'.format(rxn.id))

        try:
            rev_flux += np.exp(sol.raw.loc[rev_id])
        except:
            print('Reverse flux does not exist for {}'.format(rxn.id))

        dict_fluxes_from_sol[rxn.id] = fwd_flux - rev_flux
    fluxes_sol = pd.Series(dict_fluxes_from_sol)
    fluxes_sol = fluxes_sol.loc[reference_fluxes.index]
    ratios = fluxes_sol/(reference_fluxes/flux_scaling_ratio)

    return fluxes_sol, ratios
