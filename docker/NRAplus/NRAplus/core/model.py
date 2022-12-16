from pytfa.core.model import LCSBModel
from pytfa.io.viz import get_reaction_data
from pytfa.optim.variables import ForwardUseVariable, BackwardUseVariable

from collections import defaultdict

from cobra import Model, DictList
from cobra.core.solution import Solution

from optlang.exceptions import SolverError

from skimpy.utils.general import sanitize_cobra_vars

import pandas as pd
import numpy as np
from numpy import empty

import sympy

from NRAplus.core.solution import get_nra_solution
from NRAplus.utils.split import get_reversible_fluxes
from .variables import *
from .constraints import *

#TODO: NO TENSOR INPUT FOR CONTROL COEFFICIENTS (DATAFRAME).
#TODO: NAMESPACE OF TMODEL AND KMODEL ARE NOT THE SAME BECAUSE glc-D_e is going to be glc__D_e.
# There is a function in skimpy.utils.general - sanitize_cobra_vars which maps tmodel space against the kmodel namespace

NET = 'net'
SPLIT = 'split'

BIGM = 1000
BIGM_THERMO = 1000
MAX_REGULATION = 100  # Max allowable fold change in enzyme regulation
MAX_FOLD_CHANGE = 100

class NRAmodel(Model, LCSBModel):
    """
    NRA Model class

    NRA models allow to use Metabolic control analysis to create strain designs with
    combinatorial enzyme manipulations using constraint based analysis.
    """
    def __init__(self,
                 tmodel,
                 reference_solution=None,
                 flux_control_coefficients=None,
                 concentration_control_coefficients=None,
                 type=NET,
                 displacements=None):
        """

        :param tmodel: pyTFA model pytfa.tmodel
        :param reference_solution: Optimization solution from pyTFA Model
        :param flux_control_coefficients: flux control coeffcients from a SKiMpy models
        :param concentration_control_coefficients: concentration control coeffcients from a SKiMpy models
        :param type: NRAplus.utils.NET (default) or NRAplus.utils.SPLIT
                     to indicate the use of net or split flux control model
        """

        name = tmodel.name + '_NRA'
        LCSBModel.__init__(self, tmodel, name)

        self._cons_dict = tmodel._cons_dict.copy()
        self._var_dict = tmodel._var_dict.copy()


        self.num_enzymes = 0
        self.concentration_control_coefficients = [] #todo: what is a good way to initialize these as empty
        self.flux_control_coefficients = []
        self.reference_solution = []

        self._max_fold_enzyme_change = MAX_REGULATION
        self._max_fold_concentration_change = MAX_FOLD_CHANGE

        if type == SPLIT or type == NET:
            self.type = type
        else:
            raise ValueError("type = {} is not supported".format(type))

        if type == SPLIT and displacements is None:
            raise ValueError("To run {} analysis you need to provide pre-calculated displacements".format(type))

        self.displacements = displacements

        # self.regenerate_variables()
        if reference_solution is not None:
            # Add all the NRA constraints and variables using the reference concentrations etc.
            self.convert_model(type=type,
                               reference_solution=reference_solution,
                               flux_control_coefficients=flux_control_coefficients,
                               concentration_control_coefficients=concentration_control_coefficients,
                               displacements=displacements)

        # Repair the model dicts
        self.regenerate_variables()
        self.regenerate_constraints()
        # self.repair()

    @property
    def max_enzyme_modifications(self):
        return self.constraints['MIC_ALL'].ub + self.num_enzymes

    @max_enzyme_modifications.setter
    def max_enzyme_modifications(self, value):
        self.constraints['MIC_ALL'].ub = value - self.num_enzymes

    @property
    def max_fold_enzyme_change(self):
        return self._max_fold_enzyme_change

    @max_fold_enzyme_change.setter
    def max_fold_enzyme_change(self, value):
        """
        Set max fold change for all the enzyme regulation variables
        :param value:
        :return:
        """

        up_regulation_vars = self.get_variables_of_type(EnzymeUpRegulationVariable)
        down_regulation_vars = self.get_variables_of_type(EnzymeDownRegulationVariable)

        for var in up_regulation_vars:
            var.variable.ub = value
        for var in down_regulation_vars:
            var.variable.ub = value

        self._max_fold_enzyme_change = value

    @property
    def max_fold_concentration_change(self):
        return self._max_fold_concentration_change

    @max_fold_enzyme_change.setter
    def max_fold_concentration_change(self, value):
        """
        Set max fold change for all the enzyme regulation variables
        :param value:
        :return:
        """

        concentration_fold_vars = self.get_variables_of_type(LogConcentrationRatio)
        if np.isscalar(value):
            for var in concentration_fold_vars:
                var.variable.ub = value
                var.variable.lb = -value

        elif len(value) == 2:
            lb,ub = value
            for var in concentration_fold_vars:
                var.variable.ub = ub
                var.variable.lb = lb

        else:
            raise ValueError()

        self._max_fold_concentration_change = value

    def regenerate_variables(self):
        """
        Generates references to the cobra_model's constraints in self._var_dict
        as tab-searchable attributes of the thermo cobra_model
        :return:
        """

        _var_kinds = defaultdict(DictList)
        for k, v in self._var_dict.items():
            _var_kinds[v.__class__.__name__].append(v)

        # Add new attribute if it doesn't already exist - this ensures that the old delta_g attribute isn't overridden
        for k in _var_kinds:
            attrname = camel2underscores(k)
            if hasattr(self,attrname):
                continue
            else:
                setattr(self, attrname, _var_kinds[k])

        self._var_kinds = _var_kinds

        # Remove those attributes that are no longer in the _var_dict
        for k in self._var_kinds:
            if k in _var_kinds:
                continue
            else:
                attrname = camel2underscores(k)
                try:
                    delattr(self, attrname)
                except AttributeError:
                    pass # The attribute may not have been set up yet

    def _prune_model(self):
        """
        This hidden function removes all unnecessary variables and constraints from the scaffold thermo model
        :return:
        """
        # TODO: Please make this neater - should not have to list all types of constraints and variables, maybe choose what to keep
        variables_to_remove = ['forward_use_variable', 'backward_use_variable', 'thermo_displacement', 'min_f_lux_variable']

        for var_name in variables_to_remove:
            # Remove existing binary variables from COBRA model
            if hasattr(self, var_name):
                for var in getattr(self,var_name):
                    self._var_dict.pop(var.name)
                    self.remove_cons_vars(var.variable)
                delattr(self, var_name)

        # TODO: sort out the issue with the simultaneous use attribute being still there nra_model.simultaneous_use
        # this creates an issue when calling vars(nra_model)

        constraints_to_remove = ['forward_direction_coupling', 'backward_direction_coupling', 'forward_delta_g_coupling',
                                 'backward_delta_g_coupling', 'displacement_coupling', 'simultaneous_use', 'min_f_lux' ,
                                 'model_constraint', 'reaction_constraint']
        # Remove coupling constraints from the model
        for cons_name in constraints_to_remove:

            if hasattr(self, cons_name):
                for cons in getattr(self,cons_name):
                    self._cons_dict.pop(cons.name)
                    self.remove_cons_vars(cons.constraint)
                delattr(self,cons_name)

        # Remove default reaction flux variables from COBRA model
        for rxn in self.reactions:
            forward = rxn.forward_variable
            reverse = rxn.reverse_variable
            self.remove_cons_vars([forward, reverse])

        # Remove the mass balance constraints
        for met in self.metabolites:
            try:
                self.remove_cons_vars(self.constraints[met.id])
            except:
                print("metabolite {} does not have an associated mass balance constraint".format(met.id))

    def convert_model(self, reference_solution, flux_control_coefficients, concentration_control_coefficients, type=NET,
                      displacements=None):
        """
        This method converts the base tmodel that is used to create the NRA model into a full fledged NRA model by
        adding the new NRA variables/constraints and removing the redundant TFA constraints

        :param reference_solution: This is the output of doing tmodel.optimize() - a pandas dataframe
        :param flux_control_coefficients: This is the output of the MCA done using skimpy.
                This must be a single slice - a pandas dataframe
        :param concentration_control_coefficients: Same as the above
        :param type: NET would use NET flux control coefficients while SPLIT would use the FWD and REV fccs
        :param displacements: to use thermodynamic displacement data (for the SPLIT case)
        :return:
        """
        self.concentration_control_coefficients = concentration_control_coefficients
        self.flux_control_coefficients = flux_control_coefficients
        self.reference_solution = reference_solution

        # Split reference solution into concentrations and fluxes
        reference_concentrations_dict = {c.id: np.exp(reference_solution.loc[c.name]) for c in self.log_concentration}
        reference_concentrations = pd.Series(reference_concentrations_dict)
        reference_fluxes = get_reaction_data(self, reference_solution)

        # Remove TFA variables and constraints that are not needed
        # NOTE: this needs to be done after we use the get_reaction_data function since it uses the
        # forward and backward reaction variables which are removed when pruning
        self._prune_model()  # TODO: we still need to remove things like self.simultaneous_use

        # Add constraints/variables for NRA
        enzymes_to_regulate = concentration_control_coefficients.columns
        self._prepare_enzymes(enzymes_to_regulate)
        self._prepare_metabolites(reference_concentrations,
                                  self.concentration_control_coefficients)
        self._prepare_reactions(reference_fluxes,
                                self.flux_control_coefficients,
                                type=type,
                                displacements=displacements)
        self._prepare_thermodynamics(reference_fluxes,
                                     type=type)



    def get_solution(self):
        """
        Overrides the LCSBModel.get_solution method log fluxes

        *   :code:`solution.raw` is a clear copy of the solver output. From there one
            can access the value at solution for all the variables of the problem.
            However, looking for a reaction ID in there will only give the
            _forward_ flux. This should be used for any other variable than fluxes.

        *   :code:`solution.values` yields variables multiplied by their scaling factor
            (1 by default). Useful if you operated scaling on your equations for
            numerical reasons. This does _not_ include fluxes

        :return:
        """
        objective_value = self.solver.objective.value
        status = self.solver.status
        variables = pd.Series(data=self.solver.primal_values)

        fluxes = empty(len(self.reactions))
        rxn_index = list()
        var_primals = self.solver.primal_values

        for (i, rxn) in enumerate(self.reactions):
            rxn_index.append(rxn.id)
            # TODO access the flux vairables properly
            try:
                if self.type == NET:
                    fluxes[i] = np.exp(var_primals["LF_"+rxn.id])
                if self.type == SPLIT:
                    fluxes[i] = np.exp(var_primals["FLF_"+rxn.id])
                    if "RLF_"+rxn.id in var_primals:
                        fluxes[i] -= np.exp(var_primals["RLF_"+rxn.id])
            except KeyError:
                pass

        fluxes = pd.Series(index=rxn_index, data=fluxes, name="fluxes")

        solution = Solution(objective_value=objective_value, status=status,
                            fluxes=fluxes)

        self.solution = solution

        self.solution.raw = variables

        self.solution.values = pd.DataFrame.from_dict({k: v.unscaled
                                                       for k, v in self._var_dict.items()},
                                                       orient='index')

        return solution

    def optimize(self, objective_sense=None, **kwargs):
        """
        Call the Model.optimize function (which really is but an interface to the
        solver's. Catches SolverError in the case of no solutions. Passes down
        supplementary keyword arguments (see cobra.thermo.Model.optimize)
        :type objective_sense: 'min' or 'max'
        """

        if objective_sense:
            self.objective.direction = objective_sense
        try:
            self._optimize(self, **kwargs)
            solution = self.get_solution()
            self.solution = solution
            return solution
        except SolverError as SE:
            status = self.solver.status
            self.logger.error(SE)
            self.logger.warning('Solver status: {}'.format(status))
            raise (SE)


    """
    Hidden Utils functions
    """
    def _prepare_enzymes(self, enzymes_to_regulate):
        import numpy as np
        """
        This method adds all the new variables and constraints for the enzymes
        :return:
        """

        interventions_expression = 0
        self.num_enzymes = 0

        for rxn in self.reactions:
            if 'vmax_forward_'+rxn.id in enzymes_to_regulate:
                EUR = self.add_variable(EnzymeUpRegulationVariable, rxn, lb=0, ub=MAX_REGULATION, queue=False)
                EDR = self.add_variable(EnzymeDownRegulationVariable, rxn, lb=0, ub=MAX_REGULATION, queue=False)
                EUU = self.add_variable(EnzymeUpUseVariable, rxn, lb=0, ub=1, queue=False)
                EDU = self.add_variable(EnzymeDownUseVariable, rxn, lb=0, ub=1, queue=False)
                EUDU = self.add_variable(EnzymeUpDownUseVariable, rxn, lb=0, ub=1, queue=False)

                # Simultaneous Enzyme Regulation constraint that ensures only one of the up/down binary variables is toggled
                expression = EUU + EDU
                self.add_constraint(SimultaneousEnzymeRegulation, rxn, expression, lb=0, ub=1, queue=False)

                # Enzyme Coupling constraints : couples the binary variables to the actual up/down regulations
                expression = EUR - BIGM * EUU
                self.add_constraint(EnzymeUpCoupling, rxn, expression, lb=-BIGM, ub=-1e-09, queue=False) # TODO: Get ub from solver accuracy
                expression = EDR - BIGM * EDU
                self.add_constraint(EnzymeDownCoupling, rxn, expression, lb=-BIGM, ub=-1e-09, queue=False)

                # Intervention Coupling : denotes an intervention if either up/down regulation is toggled
                expression = EUU + EDU + BIGM * EUDU
                self.add_constraint(InterventionCoupling, rxn, expression, lb=0, ub=BIGM, queue=False)

                interventions_expression += (1 - EUDU)
                self.num_enzymes += 1

        #TODO : Might want to remove this constraint
        self.add_constraint(MaxInterventionsConstraint, hook=self, expr=interventions_expression,
                            id_='ALL', lb=0, ub=BIGM,)

    def _prepare_metabolites(self, reference_concentrations, concentration_control_coefficients):
        """
        This method adds all the new NRA constraints/variables attached to a metabolite. e.g. Concentration Response
        Balance constraints, LogConcentrationRatio variables etc.

        :param reference_concentrations: a pandas dataseries of reference concentrations
        :param concentration_control_coefficients: a dataframe of control coefficients
        :return:
        """

        if isinstance(reference_concentrations,dict):
            reference_concentrations = pd.Series(reference_concentrations)

        for met in self.metabolites:
            M = self.add_variable(kind=LogConcentrationRatio,
                                  hook=met,
                                  lb=-MAX_FOLD_CHANGE,
                                  ub=MAX_FOLD_CHANGE,
                                  queue=False)
            sanitized_id = sanitize_cobra_vars(met.id)

            if sanitized_id in concentration_control_coefficients.index:
                # Create the Concentration Response Balance constraint for the given metabolite
                control_coefficients = concentration_control_coefficients.loc[sanitized_id]
                self._create_response_balance_constraint(var=met,
                                                         nra_var=M,
                                                         control_coefficients=control_coefficients,
                                                         nra_constraint_type = ConcentrationResponseBalance)

                # Couple the concentration ratio to log concentration: ln(x_i) = LCR_i + ln(x_ref,i)
                lc = self.log_concentration.get_by_id(met.id).variable
                expression = lc - M - np.log(reference_concentrations.loc[met.id])
                self.add_constraint(kind=ConcentrationVariableCoupling,
                                    hook=met,
                                    expr=expression,
                                    lb=0,
                                    ub=0)
            else:
                print("Metabolite {} does not have an associated control coefficient".format(met.id))

    def _prepare_reactions(self, reference_fluxes, flux_control_coefficients, type=NET, displacements=None):
        """
        This method adds all the new NRA constraints/variables attached to a reaction. e.g. Flux Response
        Balance constraints, LogFluxRatio variables etc.

        :param reference_fluxes: a pandas dataseries of reference fluxes
        :param flux_control_coefficients: a pd dataframe of the flux control coefficients
        :param type : NET or SPLIT for net or fwd/bwd control coefficients
        :param displacements: Thermo displacement data for the SPLIT type to determine FWD/BWD fluxes
        :return:
        """
        # TODO: use ensure the id matiching of TFA and SKIMPY reactions names !!!!
        # implement SKIMPY_ID = sanitize_cobra_vars(TFA_ID)

        # Some book keeping
        self._NET_LOG_FLUX = dict()
        self._FWD_LOG_FLUX = dict()
        self._RVS_LOG_FLUX = dict()

        if isinstance(reference_fluxes,dict):
            reference_fluxes = pd.Series(reference_fluxes)

        # TODO SPLIT IN TWO FUNCTIONS
        for rxn in self.reactions:

            if type == NET:
                """
                Create net control-coeffcient problem
                """

                if rxn.id in flux_control_coefficients.index:

                    LFR = self.add_variable(kind=LogFluxRatio,
                                            hook=rxn,
                                            lb=-MAX_FOLD_CHANGE,
                                            ub=MAX_FOLD_CHANGE,
                                            queue=False)
                    LF  = self.add_variable(kind=LogFlux,
                                            hook=rxn,
                                            lb=-MAX_FOLD_CHANGE,
                                            ub=MAX_FOLD_CHANGE,
                                            queue=False)

                    self._NET_LOG_FLUX[rxn.id] = LF

                    control_coefficients = flux_control_coefficients.loc[rxn.id]
                    self._create_response_balance_constraint(var=rxn,
                                                             nra_var=LFR,
                                                             control_coefficients=control_coefficients,
                                                             nra_constraint_type=FluxResponseBalance)

                    # Couple the flux ratio to log flux: ln(v_i) = LFR_i + ln(x_ref,i)
                    if reference_fluxes.loc[rxn.id] > 0:
                        expression = LF - LFR - np.log(reference_fluxes.loc[rxn.id])
                    else:
                        expression = LF - LFR - np.log( -reference_fluxes.loc[rxn.id])

                    self.add_constraint(kind=FluxVariableCoupling,
                                        hook=rxn,
                                        expr=expression,
                                        lb=0,
                                        ub=0)

                else:
                    print("Flux {} does not have an associated control coefficient".format(rxn.id))

            elif type == SPLIT:
                """
                Create split fwd / bwd control-coefficient problem
                """

                forward_fluxes, reverse_fluxes = get_reversible_fluxes(reference_fluxes, displacements)

                if 'fwd_'+rxn.id in flux_control_coefficients.index:

                    FLFR = self.add_variable(ForwardLogFluxRatio, rxn, lb=-MAX_FOLD_CHANGE, ub=MAX_FOLD_CHANGE, queue=False)
                    FLF  = self.add_variable(ForwardLogFlux, rxn, lb=-MAX_FOLD_CHANGE, ub=MAX_FOLD_CHANGE, queue=False)

                    self._FWD_LOG_FLUX[rxn.id] = FLF

                    control_coefficients = flux_control_coefficients.loc['fwd_'+rxn.id]
                    self._create_response_balance_constraint(rxn, FLFR, control_coefficients, ForwardFluxResponseBalance)

                    # Couple the flux ratio to log flux: ln(v_i) = LFR_i + ln(x_ref,i)
                    expression = FLF - FLFR - np.log(forward_fluxes.loc[rxn.id])

                    self.add_constraint(ForwardFluxVariableCoupling, rxn, expression, lb=0, ub=0)
                else:
                    print("Flux {} does not have an associated forward control coefficient".format(rxn.id))

                if 'bwd_'+rxn.id in flux_control_coefficients.index and reverse_fluxes[rxn.id] > 0:
                    RLFR = self.add_variable(ReverseLogFluxRatio, rxn, lb=-MAX_FOLD_CHANGE, ub=MAX_FOLD_CHANGE, queue=False)
                    RLF  = self.add_variable(ReverseLogFlux, rxn, lb=-MAX_FOLD_CHANGE, ub=MAX_FOLD_CHANGE, queue=False)

                    self._RVS_LOG_FLUX[rxn.id] = RLF

                    control_coefficients = flux_control_coefficients.loc['bwd_'+rxn.id]
                    self._create_response_balance_constraint(rxn, RLFR, control_coefficients, ReverseFluxResponseBalance)

                    # Couple the flux ratio to log flux: ln(v_i) = LFR_i + ln(x_ref,i)
                    expression = RLF - RLFR - np.log(reverse_fluxes.loc[rxn.id])

                    self.add_constraint(ReverseFluxVariableCoupling, rxn, expression, lb=0, ub=0)
                else:
                    print("Flux {} does not have an associated reverse control coefficient".format(rxn.id))


    def _create_response_balance_constraint(self, var, nra_var, control_coefficients, nra_constraint_type, kind=NET):
        """

        :param var: reaction or metabolite to serve as hook
        :param nra_var: LogFluxRatio or LogConcentrationRatio to be used in the constraint
        :param control_coefficients: the control coefficients used in the constraint
        :param nra_constraint_type: either FluxResponseBalance or ConcentrationResponseBalance
        :param kind: either SPLIT or NET for FWD/BWD control coefficients or NET control coefficients
        :return:
        """

        # Create the associated Response Balance constraint
        expression = nra_var

        for enzyme in self.reactions:
            vmax = 'vmax_forward_' + enzyme.id

            if vmax in control_coefficients.index:
                enzyme_upregulation_var = self.variables['EU_' + enzyme.id]
                enzyme_downregulation_var = self.variables['ED_' + enzyme.id]
                expression -= control_coefficients.loc['vmax_forward_'+enzyme.id] * \
                              (enzyme_upregulation_var- enzyme_downregulation_var)  # TODO: Please check this!

        self.add_constraint(kind=nra_constraint_type,
                            hook=var,
                            expr=expression,
                            lb=0,
                            ub=0)

    def _prepare_thermodynamics(self, net_reference_fluxes, type=NET):
        """
        When kind = Net, this method will ensure that the sign of DG does not change
        When kind = Split, the method retains the old coupling and makes no changes
        :param kind:
        :param net_reference_fluxes: pandas series of fluxes
        :return:
        """

        if type == NET:

            rxns_with_delta_g = [var.id for var in self.delta_g]

            for rxn_id in rxns_with_delta_g:
                net_flux = net_reference_fluxes.loc[rxn_id]
                deltaG = self.delta_g.get_by_id(rxn_id).variable

                if net_flux > 0:
                    deltaG.ub = 0
                    deltaG.lb = -1000

                elif net_flux < 0:
                    deltaG.ub = 1000
                    deltaG.lb = 0

                else:
                    raise ValueError("Not a valid reference solution: Reaction {} has net flux = 0".format(rxn_id))

                # Make sure the _var_dict points to the updated variable
                self._var_dict[self.delta_g.get_by_id(rxn_id).name] = self.delta_g.get_by_id(rxn_id)

        elif type == SPLIT:

            epsilon = self.solver.configuration.tolerances.feasibility
            _, reverse_fluxes = get_reversible_fluxes(net_reference_fluxes, self.displacements)

            # Get all reversible fluxes with assgined delta G
            rev_rxns_with_delta_g = [var.id for var in self.delta_g if reverse_fluxes[var.id] > 0]
            irrev_rxns_with_delta_g = [var.id for var in self.delta_g if reverse_fluxes[var.id] == 0]

            for rxn_id in irrev_rxns_with_delta_g:
                net_flux = net_reference_fluxes.loc[rxn_id]
                deltaG = self.delta_g.get_by_id(rxn_id).variable

                if net_flux > 0:
                    deltaG.ub = 0
                    deltaG.lb = -1000

                elif net_flux < 0:
                    deltaG.ub = 1000
                    deltaG.lb = 0

                # Make sure the _var_dict points to the updated variable
                self._var_dict[self.delta_g.get_by_id(rxn_id).name] = self.delta_g.get_by_id(rxn_id)

            # build constraints to implement
            # if ln_v_f > ln_v_r then delta_g < 0
            # if ln_v_r < ln_v_f then delta_g > 0

            for rxn_id in rev_rxns_with_delta_g:
                FLF = self._FWD_LOG_FLUX[rxn_id]
                RLF = self._RVS_LOG_FLUX[rxn_id]

                DGR = self.delta_g.get_by_id(rxn_id).variable

                rxn = self.reactions.get_by_id(rxn_id)

                FU_rxn = self.add_variable(ForwardUseVariable, rxn, queue=False)
                BU_rxn = self.add_variable(BackwardUseVariable, rxn, queue=False)

                # create the prevent simultaneous use constraints
                # SU_rxn: FU_rxn + BU_rxn <= 1
                CLHS = FU_rxn + BU_rxn
                self.add_constraint(SimultaneousUse, rxn, CLHS, ub=1)

                # create constraints that control that log_v_f - log_v_r > 0 if FU > 0
                # UF_rxn: log_v_f - log_v_r - M FU_rxn < 0
                CLHS = FLF - RLF - FU_rxn * BIGM
                self.add_constraint(ForwardDirectionCoupling, rxn, CLHS, ub=0)

                # create constraints that control that log_v_r - log_v_f > 0 if BU > 0
                # UR_rxn: R_rxn - M RU_rxn < 0
                CLHS = RLF - FLF - BU_rxn * BIGM
                self.add_constraint(BackwardDirectionCoupling, rxn, CLHS, ub=0)

                # Note DW: implementing suggestions of maria
                # FU_rxn: 1000 FU_rxn + DGR_rxn < 1000 - epsilon
                CLHS = DGR + FU_rxn * BIGM_THERMO
                self.add_constraint(ForwardDeltaGCoupling,
                                    rxn,
                                    CLHS,
                                    ub=BIGM_THERMO - epsilon*1e3)

                # BU_rxn: 1000 BU_rxn - DGR_rxn < 1000 - epsilon
                CLHS = BU_rxn * BIGM_THERMO - DGR
                self.add_constraint(BackwardDeltaGCoupling,
                                    rxn,
                                    CLHS,
                                    ub=BIGM_THERMO - epsilon*1e3)

                # Make sure the _var_dict points to the updated variable
                self._var_dict[self.delta_g.get_by_id(rxn_id).name] = self.delta_g.get_by_id(rxn_id)

        else:
            raise ValueError("type = {} is not a valid input".format(kind))

    def _optimize(self, objective_sense=None, raise_error=False):
        """
        Optimize the model using NRA
        Parameters
        ----------
        objective_sense : {None, 'maximize' 'minimize'}, optional
            Whether fluxes should be maximized or minimized. In case of None,
            the previous direction is used.
        raise_error : bool
            If true, raise an OptimizationError if solver status is not
             optimal.
        Notes
        -----
        Only the most commonly used parameters are presented here.  Additional
        parameters for cobra.solvers may be available and specified with the
        appropriate keyword argument.
        """
        original_direction = self.objective.direction
        self.objective.direction = {"maximize": "max", "minimize": "min"}.get(
            objective_sense, original_direction
        )
        self.slim_optimize()
        solution = get_nra_solution(self, raise_error=raise_error)
        self.objective.direction = original_direction
        return solution

    def copy(self):

        # raise NotImplemented()

        from NRAplus.io.dict import nra_model_from_dict, nra_model_to_dict
        from pytfa.optim.utils import copy_solver_configuration

        dictmodel = nra_model_to_dict(self)
        new_model = nra_model_from_dict(dictmodel)

        copy_solver_configuration(self, new_model)

        return new_model
