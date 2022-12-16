# -*- coding: utf-8 -*-
"""
.. module:: skimpy
   :platform: Unix, Windows
   :synopsis: Simple Kinetic Models in Python

.. moduleauthor:: SKiMPy team

[---------]

Copyright 2017 Laboratory of Computational Systems Biotechnology (LCSB),
Ecole Polytechnique Federale de Lausanne (EPFL), Switzerland

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""

from skimpy.core.modifiers import FirstOrderSmallMoleculeModifier
from skimpy.mechanisms import *
from skimpy.utils.general import sanitize_cobra_vars
from operator import itemgetter
import re

def create_reaction_from_stoich(name,
                                met_stoich_dict,
                                model_generator,
                                inhibitors=None,
                                reaction_group=None,
                                irrev=None):

    small_molecules = model_generator.small_molecules
    reactant_relations = model_generator.reactant_relations

    this_reaction_small_molecules = TabDict([])
    this_reaction_reactants = TabDict([])

    """ 
    Split into small molecules and reactants
    and omit water 
    """

    for this_met_id, met_with_stoich in met_stoich_dict.items():
        this_met = met_with_stoich.metabolite
        stoich   = met_with_stoich.stoichiometry
        is_small_molecule = this_met.id in small_molecules
        if not (this_met.formula == WATER_FORMULA) \
           and this_met.id not in model_generator.reactants_to_exclude \
           and is_small_molecule:

            this_reaction_small_molecules[this_met_id] = stoich

        elif not (this_met.formula == WATER_FORMULA) \
           and this_met.id not in model_generator.reactants_to_exclude \
             and not is_small_molecule:

            this_reaction_reactants[this_met_id] = stoich

    # Filter inhibitors
    if inhibitors is not None:
        inhibitors = [sanitize_cobra_vars(inh.id)
                      for inh in inhibitors \
                        if not (inh.formula == WATER_FORMULA) \
                        and inh.id not in model_generator.reactants_to_exclude \
                        and not is_small_molecule]
    if inhibitors == []:
        inhibitors = None

    # TODO this currently catches transport of small molecules
    # create proper transports
    if not this_reaction_reactants:
        this_reaction_reactants = this_reaction_small_molecules
        this_reaction_small_molecules = {}

    # If there are no valid reactants
    # do not consider the reactions
    if not this_reaction_reactants:
        return None


    # Determine mechanism
    mechanism = guess_mechanism(this_reaction_reactants, inhibitors, irrev=irrev)
    # Determine reactant order
    reactants = make_reactant_set(mechanism,
                                  this_reaction_reactants,
                                  reactant_relations)
    # TODO: Can we generalize this type o?????
    #Create an inhibitor set if there are any
    if inhibitors is not None:
        inhibitor_key_value =  {'inhibitor{}'.format(i+1): e
                                for i, e in enumerate(inhibitors)}
        inhibitors = mechanism.Inhibitors(**inhibitor_key_value)

    skimpy_reaction = Reaction(name=name,
                                reactants=reactants,
                                mechanism=mechanism,
                                inhibitors=inhibitors,
                                enzyme=reaction_group,
                               )

    # Add small molecules modifiers - code changed to have 1st order modifier for IrrevMichMenten
    if 'IrrevMichaelisMenten' in mechanism.__name__:
        small_molecule_modifier = FirstOrderSmallMoleculeModifier
    else:
        small_molecule_modifier = model_generator.small_molecule_modifier

    for this_sm, stoich in this_reaction_small_molecules.items():
        small_molecule_mod = small_molecule_modifier(this_sm, stoich)
        skimpy_reaction.modifiers[small_molecule_mod.name] = small_molecule_mod

    return skimpy_reaction


def create_reaction_from_data(name,
                              reaction_data):
    """

    :param reaction_data:
    :return:
    """
    # TODO determine a data format and implement a parser to save reaction data
    # TODO E.G: GLYC: RandBiBiRevMM, S1 = A S2 = B, P1 = C , P2 = D ......
    raise NotImplmentendError


def guess_mechanism(reactants,inhibitors,irrev=None):

    stoich = [i for i in reactants.values()]
    stoich.sort(reverse=True)

    is_integer = all([float(i).is_integer() for i in stoich] )

    # Piece of code for determining if there is a transport from p to e
    is_p_e_trans = False
    compartments = [r[-1] for r in reactants.keys() if r[-1] == 'p' or r[-1] == 'e']
    unique_compartments = list(set(compartments))
    if len(unique_compartments) > 1:
        is_p_e_trans = True

    if irrev or not is_integer:
        return make_irrev_m_n_michaelis_menten(stoich)

    if inhibitors is not None:
        stoich_inhibitors = [1 for i in inhibitors]
        return make_convenience_with_inhibition(stoich,stoich_inhibitors)

    # Note this kinetic leads to a model with more parameters we will create the draft based
    # On rev hill for bibi reactions.
    # # 2) Reversible Michaelis Menten kinetics or Rev Mass Actino for p -> e transport for faster dynamics
    if stoich == [1, -1]:
        if is_p_e_trans:
            return make_rev_massaction(stoich)
        else:
            return ReversibleMichaelisMenten
    # mechanism = check_rev_michaelis_menten(reactants)
    # if mechanism is not None:
    #     return mechanism

    # 3) Test for Generalized Reversible hill
    # Stoichiometry needs to be equal to 1/-1 and same number of products
    # and substrates
    abs_eqal = lambda x: abs(x) == 1
    if all(map(abs_eqal,stoich)) and sum(stoich) == 0:
        # We use the generlized hill mechanism with h = 1 
        return make_generalized_reversible_hill_n_n_h1(stoich)
    # 4) Check for uni-bi or bi-uni
    elif stoich == [1, 1, -1]:
        return UniBiReversibleHill
    elif stoich == [1, -1, -1]:
        return BiUniReversibleHill
    # Else Conv-kinetics
    else:
        return make_convenience(stoich)



def check_boundary_reaction(cobra_reaction):
    stoich = [i for i in cobra_reaction.metabolites.values()]
    if all([i < 0 for i in stoich] ) or all([i > 0for i in stoich]):
        return True
    else:
        return False


def check_rev_michaelis_menten(reactants):
    stoich = [i for i in reactants.values()]
    stoich.sort()

    if stoich == [-1,-1,1,1]:
        return RandBiBiReversibleMichaelisMenten

    if stoich == [-1,1]:
        return ReversibleMichaelisMenten

    return None

def align_cofactors(reactants):
    """
    This function makes sure that the cofactor pairs are aligned
    :param reactants: a TabDict of reactant names as strings with the corresponding stoichiometry
    :return:
    """

    products = [k for k, v in reactants.items() if v > 0]
    substrates = [k for k, v in reactants.items() if v < 0]
    cofactor_pairs = {'nadh': 'nad',
                      'nad': 'nadh',
                      'nadp': 'nadph',
                      'nadph': 'nadp',
                      'atp': 'adp',
                      'adp': 'atp',
                      }

    sorted_products = []
    sorted_substrates = []

    # Loop through all the products
    for this_prod in products:
        prod_name = this_prod[:-2]              # the last 2 characters denote the compartment

        # Check if the product is a cofactor --> pop it from the products list --> add it to the sorted products
        if prod_name in cofactor_pairs.keys():
            products.remove(this_prod)
            sorted_products.append(this_prod)

            # Look for the other partner in substrates --> pop it from substrates --> add it to sorted substrates
            for this_substrate in substrates:
                if this_substrate[:-2] == cofactor_pairs[prod_name]:
                    substrates.remove(this_substrate)
                    sorted_substrates.append(this_substrate)

    # Add the remaining products and substrates to the sorted lists
    sorted_products.extend(products)
    sorted_substrates.extend(substrates)

    reactants_ = sorted_products + sorted_substrates
    reactants_aligned = TabDict({r:reactants[r] for r in reactants_})

    return reactants_aligned

def make_reactant_set(mechanism,
                      reactants,
                      reactant_relations):

    if 'GeneralizedReversible' in str(mechanism): # TODO please use a better way to check for Gen Rev Hill
        reactants_sorted = align_cofactors(reactants)
    else:
        reactants_sorted = TabDict(sorted(reactants.items(), key=itemgetter(1), reverse=True))

    reactant_list_temp = sorted(mechanism.reactant_stoichiometry.items(),
                           key=lambda x: (x[1],x[0]), reverse=True)

    reactant_list = [v for v,k in reactant_list_temp]

    # ToDo Sort the reactants to match the reactant list
    # using the reactant_relations

    reactants_resorted = {}
    (mets, stoichiometries) = zip(*reactants_sorted.items())
    reactants_sorted = zip(mets,stoichiometries,reactant_list)

    for met, stoich, react in reactants_sorted:
        reactants_resorted[react] = met

    return mechanism.Reactants(**reactants_resorted)

