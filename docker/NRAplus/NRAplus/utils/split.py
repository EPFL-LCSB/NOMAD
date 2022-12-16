#
def get_reversible_fluxes(net_fluxes, displacements):

    forward_fluxes = net_fluxes.copy()
    backward_fluxes = net_fluxes.copy()

    for r in net_fluxes.index:
        # Check if reaction is reversible
        if r in displacements.keys():
            forward_fluxes[r] = 1/(1-displacements[r]) * net_fluxes[r]
            backward_fluxes[r] = displacements[r]/(1-displacements[r]) * net_fluxes[r]
        else:
            forward_fluxes[r] = net_fluxes[r]
            backward_fluxes[r] = 0

    return forward_fluxes, backward_fluxes