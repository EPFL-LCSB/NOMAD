from pytfa.optim.variables import *
from pytfa.optim.constraints import *
'''
Create new constraints for NRA
'''

class FluxResponseBalance(ReactionConstraint):
    """
    Class to represent Response Balance constraints for log flux ratios.
    F - SUM (FCC_i * E) = 0
    --> ln(v/v_ref) = sum[ dlnv/dln(E/E_ref)) * ln(E/E_ref) ]
    """
    def __init__(self, reaction, expr, **kwargs):
        ReactionConstraint.__init__(self, reaction, expr, **kwargs)
    prefix = 'FRB_'

class ForwardFluxResponseBalance(ReactionConstraint):
    """
    Class to represent Response Balance constraints for Forward log flux ratios for use in the SPLIT case.
    F_Fwd - SUM (FCC_i * E) = 0
    --> ln(v_fwd/v_fwd,ref) = sum[ dlnv/dln(E/E_ref)) * ln(E/E_ref) ]
    """
    def __init__(self, reaction, expr, **kwargs):
        ReactionConstraint.__init__(self, reaction, expr, **kwargs)
    prefix = 'FFRB_'

class ReverseFluxResponseBalance(ReactionConstraint):
    """
    Class to represent Response Balance constraints for Reverse log flux ratios for use in the SPLIT case.
    F_Rev - SUM (FCC_i * E) = 0
    --> ln(v_rev/v_rev,ref) = sum[ dlnv/dln(E/E_ref)) * ln(E/E_ref) ]
    """
    def __init__(self, reaction, expr, **kwargs):
        ReactionConstraint.__init__(self, reaction, expr, **kwargs)
    prefix = 'RFRB_'

class ConcentrationResponseBalance(MetaboliteConstraint):
    """
    Class to represent Response Balance constraints for log concentration ratios.
    M - SUM (CCC_i * E) = 0
    --> ln(x/x_ref) = sum[ dlnx/dln(E/E_ref)) * ln(E/E_ref) ]
    """
    def __init__(self, metabolite, expr, **kwargs):
        MetaboliteConstraint.__init__(self, metabolite, expr, **kwargs)
    prefix = 'CRB_'

class SimultaneousEnzymeRegulation(ReactionConstraint):
    """
    Constraint that prohibits the simultaneous up and down regulation of an enzyme
    EUU + EDU <=1
    """
    def __init__(self, reaction, expr, **kwargs):
        ReactionConstraint.__init__(self, reaction, expr, **kwargs)
    prefix = 'SER_'

class EnzymeUpCoupling(ReactionConstraint):
    """
    Constraint to couple the up-regulation binary variable to the value of the enzyme up-regulation
    EU - BigM*EUU <=0
    """
    def __init__(self, reaction, expr, **kwargs):
        ReactionConstraint.__init__(self, reaction, expr, **kwargs)
    prefix = 'EUC_'

class EnzymeDownCoupling(ReactionConstraint):
    """
    Constraint to couple the down-regulation binary variable to the value of the enzyme down-regulation
    ED - BigM*EDU <=0
    """
    def __init__(self, reaction, expr, **kwargs):
        ReactionConstraint.__init__(self, reaction, expr, **kwargs)
    prefix = 'EDC_'

class EnzymeDeltaG(ReactionConstraint):
    """
    Class to represent thermodynamics constraints.

    G - DGoRerr_Rxn + RT * sum (M-tilda)   = 0
    M-tilda = M + ln(x_ref)
    """

    prefix = 'EG_'

class InterventionCoupling(ReactionConstraint):
    """
    Class to couple the up/down regulation with the number of interventions
    EU + ED + BIGM*EUDU <= BIGM

    i.e. if there is up/down regulation of the enzyme, this counts as an intervention
    """
    def __init__(self, reaction, expr, **kwargs):
        ReactionConstraint.__init__(self, reaction, expr, **kwargs)
    prefix = 'IC_'

class MaxInterventionsConstraint(GenericConstraint):
    """
        Class to limit the number of interventions
        Sum of (1 - EUDU) < a given number

    """
    def __init__(self, model, expr, id_, **kwargs):

        GenericConstraint.__init__(self,
                                   id_=id_,
                                   expr=expr,
                                   model=model,
                                   **kwargs)
    prefix = 'MIC_'

class ConcentrationVariableCoupling(MetaboliteConstraint):
    """
    Class to represent the coupling of the Log concentration to the LogConcentrationRatio and the reference concentration
    ln(x) = M + ln(x_ref)
    """
    def __init__(self, metabolite, expr, **kwargs):
        MetaboliteConstraint.__init__(self, metabolite, expr, **kwargs)
    prefix = 'CVC_'

class FluxVariableCoupling(ReactionConstraint):
    """
    Class to represent the coupling of the Log concentration to the LogConcentrationRatio and the reference concentration
    ln(v) = F + ln(v_ref)
    """
    def __init__(self, reaction, expr, **kwargs):
        ReactionConstraint.__init__(self, reaction, expr, **kwargs)
    prefix = 'FVC_'

class ForwardFluxVariableCoupling(ReactionConstraint):
    """
    Class to represent the coupling of the Log concentration to the LogConcentrationRatio and the reference concentration
    ln(v) = F + ln(v_ref)
    """
    def __init__(self, reaction, expr, **kwargs):
        ReactionConstraint.__init__(self, reaction, expr, **kwargs)
    prefix = 'FFVC_'


class ReverseFluxVariableCoupling(ReactionConstraint):
    """
    Class to represent the coupling of the Log concentration to the LogConcentrationRatio and the reference concentration
    ln(v) = F + ln(v_ref)
    """
    def __init__(self, reaction, expr, **kwargs):
        ReactionConstraint.__init__(self, reaction, expr, **kwargs)
    prefix = 'RFVC_'

class ObjectiveConstraint(GenericConstraint):
    """
    Class to represent a constraint on the objective function value. This is valuable during variability analysis
    especially when the objective is an expression.
    eg. OBJ = LFR_TYRTA - LFR_GLCtex
    max obj = 0.15
    0.1 <= LFR_TYRTA-LFR_GLCtex <= 0.2
    """

    def __init__(self, model, expr, id_, **kwargs):

        GenericConstraint.__init__(self,
                                   id_=id_,
                                   expr=expr,
                                   model=model,
                                   **kwargs)

    prefix = 'OC_'


