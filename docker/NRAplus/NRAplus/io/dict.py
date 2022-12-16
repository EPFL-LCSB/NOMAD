# -*- coding: utf-8 -*-
from __future__ import absolute_import

from collections import OrderedDict
from operator import attrgetter, itemgetter

from numpy import bool_, float_
from six import iteritems, string_types

from cobra.core import Gene, Metabolite, Model, Reaction
from cobra.util.solver import set_objective

from pytfa.optim.variables import ReactionVariable, MetaboliteVariable, ModelVariable
from pytfa.optim.constraints import ReactionConstraint, MetaboliteConstraint, ModelConstraint, GenericConstraint

from NRAplus.core.model import NRAmodel

from optlang.util import expr_to_json, parse_expr

# TODO: These are all just modifications of the cobra code, do I need to write separate functions?
#  Or just overwrite the global variables?
_REQUIRED_REACTION_ATTRIBUTES = [
    "id", "name", "metabolites", "lower_bound", "upper_bound",
    "gene_reaction_rule"]
_ORDERED_OPTIONAL_REACTION_KEYS = [
    "subsystem", "notes", "annotation"]
_OPTIONAL_REACTION_ATTRIBUTES = {
    "subsystem": "",
    "notes": {},
    "annotation": {},
}
BASE_NAME2HOOK = {
                    ReactionVariable    :'reactions',
                    ReactionConstraint  :'reactions',
                    MetaboliteVariable  :'metabolites',
                    MetaboliteConstraint:'metabolites',
                }

SOLVER_DICT = {
    'optlang.gurobi_interface':'optlang-gurobi',
    'optlang.cplex_interface':'optlang-cplex',
    'optlang.glpk_interface':'optlang-glpk',
}
_ORDERED_OPTIONAL_MODEL_KEYS = ["name", "compartments", "notes", "annotation"]
_OPTIONAL_MODEL_ATTRIBUTES = {
    "name": None,
    #  "description": None, should not actually be included
    "compartments": [],
    "notes": {},
    "annotation": {},
}

def reaction_to_dict(reaction):
    """
    Method to convert reaction information into a dictionary. Copied from the cobra function of the same name with
    a few modifications

    :param reaction:
    :return: Ordered dictionary of reactions with relevant information
    """
    from cobra.io.dict import metabolite_to_dict, gene_to_dict, _update_optional, _fix_type

    new_reaction = OrderedDict()
    for key in _REQUIRED_REACTION_ATTRIBUTES:
        if key != "metabolites":
            new_reaction[key] = _fix_type(getattr(reaction, key))
            continue
        mets = OrderedDict()
        for met in sorted(reaction.metabolites, key=attrgetter("id")):
            mets[str(met)] = reaction.metabolites[met]
        new_reaction["metabolites"] = mets
    _update_optional(reaction, new_reaction, _OPTIONAL_REACTION_ATTRIBUTES,
                     _ORDERED_OPTIONAL_REACTION_KEYS)
    return new_reaction

def nra_model_to_dict(model, sort=False):
    """
    Method that copies all the information from an NRA model into a dictionary

    :param model: NRA model to be copied
    :return: obj: Dictionary with all the stored information
    """
    from cobra.io.dict import metabolite_to_dict, gene_to_dict, _update_optional, _fix_type
    from pytfa.io.dict import get_solver_string, \
                                obj_to_dict, var_to_dict, cons_to_dict, \
                                _add_thermo_metabolite_info, _add_thermo_reaction_info

    # Take advantage of cobra's dict serialization for metabolites and reactions
    # Copied from obj = cbd.model_to_dict(model)
    obj = OrderedDict()
    obj["metabolites"] = list(map(metabolite_to_dict, model.metabolites))
    obj["reactions"] = list(map(reaction_to_dict, model.reactions))
    obj["genes"] = list(map(gene_to_dict, model.genes))
    obj["id"] = model.id
    _update_optional(model, obj, _OPTIONAL_MODEL_ATTRIBUTES,
                     _ORDERED_OPTIONAL_MODEL_KEYS)
    if sort:
        get_id = itemgetter("id")
        obj["metabolites"].sort(key=get_id)
        obj["reactions"].sort(key=get_id)
        obj["genes"].sort(key=get_id)

    obj['solver'] = get_solver_string(model)
    obj['objective'] = obj_to_dict(model)

    # Copy variables, constraints
    obj['variables'] = list(map(var_to_dict, model._var_dict.values()))
    obj['constraints'] = list(map(cons_to_dict, model._cons_dict.values()))

    # Thermo data and variables
    obj['kind'] = 'ThermoModel'
    obj['thermo_data'] = model.thermo_data  # it's a dict
    obj['name'] = model.name
    obj['temperature'] = model.TEMPERATURE
    obj['min_ph'] = model.MIN_pH
    obj['max_ph'] = model.MAX_pH
    is_thermo = True

    # Metabolite and Reaction-level cleanup
    for rxn_dict in obj['reactions']:
        rxn = model.reactions.get_by_id(rxn_dict['id'])

        if is_thermo:
            _add_thermo_reaction_info(rxn, rxn_dict)

    # Peptides and Thermo
    for met_dict in obj['metabolites']:
        the_met_id = met_dict['id']
        is_peptide = False

        if is_thermo and not is_peptide:  # peptides have no thermo
            the_met = model.metabolites.get_by_id(the_met_id)
            _add_thermo_metabolite_info(the_met, met_dict)
            met_dict['kind'] = 'Metabolite'

    # Relaxation info
    try:
        obj['relaxation'] = model.relaxation
    except AttributeError:
        pass

    # Control coefficients
    obj['concentration_control_coefficients'] = model.concentration_control_coefficients.copy()
    obj['flux_control_coefficients'] = model.flux_control_coefficients.copy()
    obj['reference_solution'] = model.reference_solution.copy()

    # Other NRA attributes
    obj["max_enzyme_modifications"] = model.max_enzyme_modifications
    obj["max_fold_enzyme_change"] = model.max_fold_enzyme_change
    obj["max_fold_concentration_change"] = model.max_fold_concentration_change
    obj["num_enzymes"] = model.num_enzymes

    return obj

def cobra_model_from_dict(obj):
    """
    Copied from the cobra.io.dict file. Major difference is when adding reactions, we don't add certain
    reaction constraints.

    :param obj: dictionary containing all the old NRA model's constraints and variables
    :return:
    """
    from cobra.io.dict import metabolite_from_dict, reaction_from_dict, gene_from_dict

    if 'reactions' not in obj:
        raise ValueError('Object has no reactions attribute. Cannot load.')
    cobra_model = Model()
    cobra_model.add_metabolites(
        [metabolite_from_dict(metabolite) for metabolite in obj['metabolites']]
    )
    cobra_model.genes.extend([gene_from_dict(gene) for gene in obj['genes']])
    cobra_model.add_reactions(
        [reaction_from_dict(reaction, cobra_model) for reaction in obj['reactions']]
    )

    for k, v in iteritems(obj):
        if k in {'id', 'name', 'notes', 'compartments', 'annotation'}:
            setattr(cobra_model, k, v)

    return cobra_model

def tfa_model_from_dict(obj, solver=None, custom_hooks=None):
    """
    Creates a skeleton tfa model to be used in creating the subsequent NRA model
    :param obj: dictionary with all the information from the NRA model to be copied
    :param solver:
    :param custom_hooks: hooks to be associated with the correct constraint/variable class
    :return: tmodel: a skeleton thermodynamic model
    """
    from pytfa.thermo.tmodel import ThermoModel
    from pytfa.io.dict import init_thermo_model_from_dict
    from pytfa.io.dict import add_custom_classes, get_model_constraint_subclasses, get_model_variable_subclasses

    # Create a cobra model from the dictionary
    cobra_model = cobra_model_from_dict(obj)

    if solver is not None:
        try:
            cobra_model.solver = solver
        except SolverNotFound as snf:
            raise snf
    else:
        try:
            cobra_model.solver = obj['solver']
        except KeyError:
            pass

    # Create a thermo model from the cobra model and dictionary
    # TODO: Should be a condition to check if it is a thermomodel?
    tmodel = ThermoModel(thermo_data=obj['thermo_data'],
                         model=cobra_model,
                         name=obj['name'],
                         temperature=obj['temperature'],
                         min_ph=obj['min_ph'],
                         max_ph=obj['max_ph'])

    # Update the thermo data for all the reactions and metabolites
    tmodel = init_thermo_model_from_dict(tmodel, obj)

    tmodel._push_queue()
    tmodel.prepare() # TODO: Is this necessary?

    tmodel.repair()

    return tmodel

def nra_model_from_dict(obj, solver=None, custom_hooks=None):
    """
    Part of the self.copy function: Creates a new NRA model from a dictoinary containing all the info from the old NRA
    model to be copied
    :param obj: dictionary containing all the information from the old NRA model to be copied
    :param solver:
    :param custom_hooks: hooks that connect to the relevant classes of variables/constraints
    :return: nra_model: the copied NRA model
    """
    from pytfa.io.dict import add_custom_classes, make_subclasses_dict, get_model_constraint_subclasses, get_model_variable_subclasses
    from pytfa.io.dict import rebuild_obj_from_dict

    if custom_hooks is None:
        custom_hooks = dict()

    custom_hooks.update(BASE_NAME2HOOK)

    # Create a thermo model from the dictionary
    tmodel = tfa_model_from_dict(obj, solver, custom_hooks)

    # Create the NRA model from the thermo model
    nra_model = NRAmodel(tmodel)
    nra_model._prune_model() # This pruning removes all the redundant TFA variables/constraints
    name2class, name2hook = add_custom_classes(nra_model, custom_hooks)

    # Update all the NRA variables
    for the_var_dict in obj['variables']:
        this_id = the_var_dict['id']
        classname = the_var_dict['kind']
        lb = the_var_dict['lb']
        ub = the_var_dict['ub']
        scaling_factor = the_var_dict['scaling_factor']
        this_var_name = the_var_dict['name']

        if  this_var_name in nra_model._var_dict:
            nra_model.variables[this_var_name].lb = lb
            nra_model.variables[this_var_name].ub = ub
            print("{} already exists".format(this_var_name))

        elif classname in get_model_variable_subclasses(): #TODO: I am still not sure where ModelConstraints will get handled
            hook = nra_model
            this_class = get_model_variable_subclasses()[classname]
            nv = nra_model.add_variable(kind=this_class,
                                        hook=hook,
                                        id_=this_id,
                                        ub=ub,
                                        lb=lb,
                                        queue=False,)

        elif classname in name2class:
            hook = name2hook[classname].get_by_id(this_id)
            this_class = name2class[classname]
            nv = nra_model.add_variable(kind=this_class,
                                        hook=hook,
                                        ub=ub,
                                        lb=lb,
                                        queue=False,)

        else:
            print('Class {} serialization not handled yet'.format(classname)) # TODO: Make this raise an error

    nra_model._push_queue() #Doesn't seem to work!?
    variable_parse_dict = {x.name: x for x in nra_model.variables}

    # Update all the NRA constraints
    for the_cons_dict in obj['constraints']:
        this_id = the_cons_dict['id']
        classname = the_cons_dict['kind']
        new_expr = parse_expr(the_cons_dict['expression'],
                              local_dict=variable_parse_dict)
        this_cons_name = the_cons_dict["name"]
        lb = the_cons_dict['lb']
        ub = the_cons_dict['ub']

        # Look for the corresponding class:
        if this_cons_name in nra_model._cons_dict:
            nra_model.constraints[this_cons_name].ub = ub
            nra_model.constraints[this_cons_name].lb = lb

        elif classname in get_model_constraint_subclasses():
            hook = nra_model
            this_class = get_model_constraint_subclasses()[classname]
            nc = nra_model.add_constraint(kind=this_class, hook=hook,
                                          expr=new_expr, id_=this_id,
                                          ub=ub,
                                          lb=lb,
                                          queue=False)

        elif classname in name2class:
            hook = name2hook[classname].get_by_id(this_id)
            this_class = name2class[classname]
            nc = nra_model.add_constraint(kind=this_class, hook=hook,
                                          expr=new_expr,
                                          ub=ub,
                                          lb=lb,
                                          queue=False)

        elif classname in make_subclasses_dict(GenericConstraint): # This handles the MIC constraint
            hook = nra_model
            this_class = make_subclasses_dict(GenericConstraint)[classname]
            nc = nra_model.add_constraint(kind=this_class, hook=hook,
                                          expr=new_expr,
                                          id_=this_id,
                                          lb=lb,
                                          ub=ub,
                                          queue=False)

        else:
            print('Class {} serialization not handled yet'.format(classname))

    nra_model.repair()

    # Get the objective function
    try:
        rebuild_obj_from_dict(nra_model, obj['objective'])
    except KeyError:
        pass

    # Relaxation info
    try:
        nra_model.relaxation = obj['relaxation']
    except KeyError:
        pass

    # Other attributes TODO: Should we loop through this?
    nra_model.num_enzymes = obj["num_enzymes"]
    nra_model.max_fold_enzyme_change = obj["max_fold_enzyme_change"]
    nra_model.max_fold_concentration_change = obj["max_fold_concentration_change"]
    nra_model.max_enzyme_modifications = obj["max_enzyme_modifications"]
    nra_model.concentration_control_coefficients = obj["concentration_control_coefficients"]
    nra_model.flux_control_coefficients = obj["flux_control_coefficients"]
    nra_model.reference_solution = obj["reference_solution"]

    nra_model.regenerate_variables()
    nra_model.regenerate_constraints()

    return nra_model
