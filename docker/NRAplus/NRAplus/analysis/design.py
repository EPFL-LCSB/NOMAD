import pandas as pd
import numpy as np

def get_modified_reactions(solution,nra_model):
    """

    :param solution: pd.Series
    :param nra_model:
    :return:
    """
    return [v.id for v in nra_model.enzyme_up_down_use_variable
            if solution[v.name] == 0]

def get_modifications(solution,nra_model,reactions):
    """

    :param solution: pd.Series
    :param nra_model:
    :param reactions:
    :return:
    """

    modifications = []
    up = nra_model.enzyme_up_regulation_variable.get_by_any(reactions)
    down = nra_model.enzyme_down_regulation_variable.get_by_any(reactions)

    epsilon = nra_model.solver.configuration.tolerances.feasibility

    for u, d in zip(up, down):
        up_regulation = solution[u.name]
        down_regulation = solution[d.name]

        if up_regulation > 0 and down_regulation < epsilon:
            regulation = 'up'
            change = up_regulation
        elif down_regulation > 0 and up_regulation < epsilon:
            regulation = 'down'
            change = -down_regulation
        else:
            regulation = 'no-change'
            change = 0

        this_mod = {
            'regulation': regulation,
            'fold_change': change,
            'change': np.exp(change)
        }
        modifications.append(this_mod)

    return modifications

def get_desing_from_solution(solution, nra_model):
    """

    :param solution: nra_solution object or pd.Series
    :param nra_model:
    :return:
    """
    #TODO Check here for input type if solution object or dataframe

    # Get the enzyme modifications
    modified_reactions = get_modified_reactions(solution.raw, nra_model)

    # Get which enzyme modification
    modifications = get_modifications(solution.raw, nra_model, modified_reactions)

    # make a table
    table = pd.DataFrame(modifications, index=modified_reactions)

    return table

def enumerate_solutions(model, max_iter = 50, threshold = 0.8):
    """
    Takes an NRAmodel and enumerates the possible combinations of enzyme changes
    :param model: The NRA model that is fed in
    :param max_iter: Maximum number of iterations allowed
    :param threshold: Min objective - a value of 0.8 means we must be within 20% of the maximal objective value

    :return: profiles: A list of raw solutions for each possible profile
    :return: list_of_designs: A list of dictionaries - each dictionary consists of a design with the keys being the
                            manipulated enzymes and the values being the fold change. In addition, there is a 'solution'
                            key that stores the NRA solution value as well
    :return: nra_model: The new nra_model with all the new constraints forbidding the extant profiles
    """
    from pytfa.optim.constraints import ForbiddenProfile
    from optlang.exceptions import SolverError
    from optlang.interface import INFEASIBLE

    nra_model = model.copy()
    solution = nra_model.optimize()
    min_obj_val = solution.objective_value * threshold

    max_regulations = nra_model.max_enzyme_modifications
    profiles = dict()
    iter_count = 0
    epsilon = nra_model.solver.configuration.tolerances.feasibility

    # # Enzyme upregulation and downregulation variables (EU + ED)
    # enzyme_downregulation_vars = [v.variable for v in nra_model.enzyme_down_use_variable]
    # enzyme_upregulation_vars = [v.variable for v in nra_model.enzyme_up_use_variable]
    #
    # # Sum of all EU + DU vars
    # expr = sum(enzyme_upregulation_vars) + sum(enzyme_downregulation_vars)
    # nra_model.add_constraint(ForbiddenProfile,
    #                            hook = nra_model,
    #                            expr = expr,
    #                            id_ = 'BRDS_active',
    #                            lb = 0,
    #                            ub = max_regulations,) #TODO: check this please
    list_of_designs = []
    while nra_model.solver.status != INFEASIBLE and iter_count < max_iter and solution.objective_value >= min_obj_val:
        try:
            solution = nra_model.optimize() #TODO: check this optimize()
            print("The solution is {}".format(solution))
        except SolverError:
            break
            print('inf')

        profiles[iter_count] = solution # Changed to the solution object itself from .raw

        if iter_count > 0:
            steady_state_error = sum((profiles[iter_count-1].raw - profiles[iter_count].raw)**2)
            print("The difference between the profiles is : {}".format(steady_state_error))
        else:
            sse =0

        nra_model.logger.debug(str(iter_count) + ' - ' + str(sse))

        # Get all enzymes that are active in this iteration
        active_enzyme_regulations = [v for v in nra_model.enzyme_up_use_variable
                                     if np.abs(nra_model.solver.primal_values[v.name] - 1) <= epsilon]
        active_enzyme_regulations.extend([v for v in nra_model.enzyme_down_use_variable
                                     if np.abs(nra_model.solver.primal_values[v.name] - 1) <= epsilon])

        dict_design = {}
        for e in active_enzyme_regulations:
            e_name = e.name.replace('U_','_')
            print("Enzyme {} regulation --> {}".format(e_name,nra_model.solver.primal_values[e_name]))
            dict_design[e_name] = nra_model.solver.primal_values[e_name]

        dict_design['solution'] = solution.objective_value
        list_of_designs.append(dict_design)

        bool_id = ''.join('-' + v.id for v in active_enzyme_regulations)

        # Make the expression to forbid this expression profile to happen again
        # FP-ENO_GLYCK2_FBP: EU/ED_ENO + EU/ED_GLYCK2 + EU/ED_FBP <= 3-1 = 2
        expr = sum(active_enzyme_regulations)

        try:
            nra_model.add_constraint(ForbiddenProfile,
                                       hook = nra_model,
                                       expr = expr,
                                       id_ = str(iter_count) + bool_id,
                                       lb = 0,
                                       ub = len(active_enzyme_regulations) - 1)

            iter_count += 1
        except:
            ValueError("Constraints can no longer be added")
            break
    return profiles, list_of_designs, nra_model

def enzyme_variability_analysis(model, sol, threshold = 0.05):
    """
    This function conducts variability analysis for the enzyme activities of a given design.

    :param model: An existing NRA model with an objective function defined already
    :param sol: an existing NRA solution
    :param threshold: The percentage of flexibility afforded to the objective function 0.05 --> 5% +/- the objective value
    :return:
            df_variability: A pandas dataframe containing min and max enzyme activity levels
            for each enzyme in a given design
    """
    from NRAplus.core.constraints import ObjectiveConstraint

    # Get the existing design
    active_enzyme_regulations = sol.raw.loc[[i for i in sol.raw.index if ('EU_' in i or 'ED_' in i)
                                             and sol.raw.loc[i] >= model.solver.configuration.tolerances.integrality]]

    # Create a copy of the NRA model and introduce a constraint that fixes the objective value
    nra_model = model.copy()
    expr = nra_model.objective.expression
    lb = (1-threshold)*sol.objective_value
    ub = (1+threshold)*sol.objective_value
    nra_model.add_constraint(ObjectiveConstraint,
                             hook=nra_model,
                             expr=expr,
                             id_='OC',
                             lb= lb,
                             ub=ub )

    # Ensure that ALL the enzymes in the design are turned on and have unlimited activity
    for this_enzyme in active_enzyme_regulations.keys():
        nra_model.variables[this_enzyme].lb = 1e-3
        nra_model.variables[this_enzyme].ub = 100

    df_variability = {}

    # Determine the min, and max values of each enzyme activity
    for this_enzyme in active_enzyme_regulations.keys():

        nra_model.objective = nra_model.variables[this_enzyme]
        nra_model.objective_direction = 'min'
        this_sol_min = nra_model.optimize()

        nra_model.objective_direction = 'max'
        this_sol_max = nra_model.optimize()

        dict_temp = {'min': this_sol_min.objective_value,
                     'max': this_sol_max.objective_value}
        df_variability[this_enzyme] = dict_temp

    df_variability = pd.DataFrame(df_variability).T

    return df_variability
