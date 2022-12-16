from cobra.core.solution import Solution
from cobra.util.solver import check_solver_status

from numpy import empty, nan
from pandas import  Series

def get_nra_solution(model, raise_error=False):
    """
    Generate a solution representation of the current solver state.
    Parameters
    ---------
    model : cobra.Model
        The model whose reactions to retrieve values for.
    raise_error : bool
        If true, raise an OptimizationError if solver status is not optimal.
    Returns
    -------
    cobra.Solution
    Note
    ----
    This is only intended for the `optlang` solver interfaces and not the
    legacy solvers.
    """
    check_solver_status(model.solver.status, raise_error=raise_error)

    rxn_index = list([rxn.id for rxn in model.reactions])
    # met_index = list([met.id for met in model.metabolites])

    reactions = model.reactions
    metabolites = model.metabolites

    fluxes = empty(len(reactions))
    # reduced = empty(len(reactions))
    # shadow = empty(len(metabolites))

    # TODO Make an NRA solution object that is more useful for the user
    return Solution(
        model.solver.objective.value,
        model.solver.status,
        Series(index=rxn_index, data=fluxes, name="fluxes"),
        # Series(index=rxn_index, data=reduced, name="reduced_costs"),
        # Series(index=met_index, data=shadow, name="shadow_prices"),
    )