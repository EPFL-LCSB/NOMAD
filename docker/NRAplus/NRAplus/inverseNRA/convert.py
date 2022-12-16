from cobra.util.solver import interface_to_str

from pytfa.optim.utils import symbol_sum
from pytfa.io.viz import get_reaction_data

from NRAplus.utils.split import get_reversible_fluxes
from NRAplus.core.variables import LogFlux, ForwardLogFlux, ReverseLogFlux, LogConcentration
from NRAplus.core.model import SPLIT, NET
from NRAplus.inverseNRA.variables import *
from NRAplus.inverseNRA.constraints import *

import numpy as np

BOTH = 'fluxes_concentrations'
FLUXES = 'fluxes'
CONCENTRATIONS = 'concentrations'

# TODO
CPLEX = 'cplex'
GUROBI = 'gurobi'

def convert_to_inverse(nra_model, tmodel, target_solution, inplace=True, type=FLUXES):
    """

    :param nra_model:
    :param tmodel:
    :param target_solution:
    :param inplace:
    :param type:
    :return:
    """

    if inplace:
        inra_model = nra_model
    else:
        inra_model = nra_model.copy()

    # Change solver configuration
    solver_str = interface_to_str(nra_model.solver.interface)
    if solver_str == CPLEX:
        pass
        #TODO: Find the parameter for cplex
    elif solver_str == GUROBI:
        inra_model.solver.problem.Params.NonConvex = 2
    else:
        raise ValueError('Solver does most likly not support non convex MIQP (If yes please implement)')

    # Add expression for the distance
    #  (log_v - log_v_target)**2 = 1.0*log_v*log_v - log_v_target*log_v + log_v_target*log_v_target

    # NOTE QUADRATIC CONSTRAINTS NEED TO BE WRITTEN IN IN THE FOLLOWING FORMAT:
    # NUMBER * VAR1 * VAR2

    # Fetch all variable that are in the solution
    distance_fluxes = 0
    if type == FLUXES or type == BOTH:
        reference_fluxes = get_reaction_data(tmodel, target_solution)

        if inra_model.type == SPLIT:
            log_forward_variables = inra_model.get_variables_of_type(ForwardLogFlux)
            log_reverse_variables = inra_model.get_variables_of_type(ReverseLogFlux)

            forward_fluxes, reverse_fluxes = get_reversible_fluxes(reference_fluxes, nra_model.displacements)

            forward_vars_and_fluxes = [(log_forward_variables.get_by_id(r.id), np.log(forward_fluxes[r.id]))
                                       for r in inra_model.forward_log_flux]

            reverse_vars_and_fluxes = [(log_reverse_variables.get_by_id(r.id), np.log(reverse_fluxes[r.id]))
                                       for r in inra_model.reverse_log_flux]

            vars_and_fluxes = forward_vars_and_fluxes + reverse_vars_and_fluxes

        elif inra_model.type == NET:
            log_flux = inra_model.get_variables_of_type(LogFlux)
            # TODO MAKE USE OF THE SIGN FUNCTION CHECK AGAIN IN MAIN CODE!!!!!
            vars_and_fluxes = [(log_flux.get_by_id(r.id), np.log(np.abs(reference_fluxes[r.id])))
                               for r in inra_model.log_flux]
            
        # QUADRATIC TERMS
        fluxes_quadratic = symbol_sum([1.0*v.variable*v.variable for v, v_ref in vars_and_fluxes])

        # LINEAR TERMS
        fluxes_linear = symbol_sum([2.0*v.variable*v_ref for v, v_ref in vars_and_fluxes])

        # CONSTANT
        fluxes_constant = symbol_sum([v_ref*v_ref for v, v_ref in vars_and_fluxes])

        # UPDATE
        distance_fluxes += fluxes_quadratic - fluxes_linear + fluxes_constant

    distance_concentrations = 0
    if type == CONCENTRATIONS or type == BOTH:

        reference_concentrations_dict = {c.id: reference_solution.loc[c.name]
                                         for c in inra_model.log_concentration}
        reference_concentrations = pd.Series(reference_concentrations_dict)

        log_concentration = inra_model.get_variables_of_type(LogConcentration)
        vars_and_concs = [(log_concentration.get_by_id(m    .id), reference_concentrations[m.id])
                          for m in inra_model.metabolites]

        # QUADRATIC TERMS
        concentrations_quadratic = symbol_sum([1.0*v.variable*v.variable for v, v_ref in vars_and_concs])

        # LINEAR TERMS
        concentrations_linear = symbol_sum([2.0*v.variable*v_ref for v, v_ref in vars_and_concs])

        # CONSTANT
        concentrations_constant = symbol_sum([v_ref*v_ref for v, v_ref in vars_and_concs])

        # UPDATE
        distance_concentrations += concentrations_quadratic - concentrations_linear + concentrations_constant

    distance = distance_concentrations + distance_fluxes

    ## TODO OPTLANG ERROR:  NotImplementedError: Quadratic constraints ( ... ) are not supported yet
    # # Add distance variable
    # D = inra_model.add_variable(TargetDistance, hook=inra_model, id_='DISTANCE')
    #
    # expression = D - distance
    # # Add distance coupling constraint
    # inra_model.add_constraint(TargetDistanceCoupling,
    #                           hook=inra_model, expr=expression,
    #                           id_='ALL', lb=0, ub=0,)

    # Set objective to distance + direction to min
    inra_model.objective = distance
    inra_model.objective_direction = 'min'

    return inra_model
