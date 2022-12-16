from pytfa.optim.variables import *
from pytfa.optim.constraints import *
'''
New Variables for NRA
'''

class LogFlux(ReactionVariable):
    """
    Class to represent the relative flux ln(v)
    """
    def __init__(self,reaction,**kwargs):
        ReactionVariable.__init__(self, reaction, **kwargs)
    prefix = 'LF_'

class LogFluxRatio(ReactionVariable):
    """
    Class to represent the relative flux ln(v/v_ref)
    """
    def __init__(self,reaction,**kwargs):
        ReactionVariable.__init__(self, reaction, **kwargs)
    prefix = 'LFR_'

class ForwardLogFlux(ReactionVariable):
    """
    Class to represent the relative flux ln(v)
    """
    def __init__(self,reaction,**kwargs):
        ReactionVariable.__init__(self, reaction, **kwargs)
    prefix = 'FLF_'

class ForwardLogFluxRatio(ReactionVariable):
    """
    Class to represent the relative flux ln(v/v_ref)
    """
    def __init__(self,reaction,**kwargs):
        ReactionVariable.__init__(self, reaction, **kwargs)
    prefix = 'FLFR_'

class ReverseLogFlux(ReactionVariable):
    """
    Class to represent the relative flux ln(v)
    """
    def __init__(self,reaction,**kwargs):
        ReactionVariable.__init__(self, reaction, **kwargs)
    prefix = 'RLF_'

class ReverseLogFluxRatio(ReactionVariable):
    """
    Class to represent the relative flux ln(v/v_ref)
    """
    def __init__(self,reaction,**kwargs):
        ReactionVariable.__init__(self, reaction, **kwargs)
    prefix = 'RLFR_'


class EnzymeUpRegulationVariable(ReactionVariable):
    """
    Class to represent the relative enzyme activity ln(E/E_ref) when E > E_ref
    """
    def __init__(self,reaction,**kwargs):
        ReactionVariable.__init__(self, reaction, **kwargs)
    prefix = 'EU_'

class EnzymeDownRegulationVariable(ReactionVariable):
    """
    Class to represent the relative enzyme activity ln(E/E_ref) where E < E_ref
    """
    def __init__(self,reaction,**kwargs):
        ReactionVariable.__init__(self, reaction, **kwargs)
    prefix = 'ED_'

class LogConcentrationRatio(MetaboliteVariable):
    """
    Class to represent the relative metabolite concentration ln(x/x_ref)
    """
    def __init__(self,metabolite,**kwargs):
        MetaboliteVariable.__init__(self, metabolite, **kwargs)
    prefix = 'LCR_'


class EnzymeUpUseVariable(ReactionVariable, BinaryVariable):
    """
    Class to represent a binary up-regulation use variable. It is used to
    activate enzyme up regulation
    """

    def __init__(self, reaction, **kwargs):
        ReactionVariable.__init__(self, reaction,
                                  type=get_binary_type(),
                                  **kwargs)

    prefix = 'EUU_'

class EnzymeDownUseVariable(ReactionVariable, BinaryVariable):
    """
    Class to represent a binary down-regulation use variable. It is used to
    activate enzyme down-regulation
    """

    def __init__(self, reaction, **kwargs):
        if not 'lb' in kwargs:
            kwargs['lb'] = 0
        if not 'ub' in kwargs:
            kwargs['ub'] = 1

        ReactionVariable.__init__(self, reaction,
                                  type=get_binary_type(),
                                  **kwargs)

    prefix = 'EDU_'

class EnzymeUpDownUseVariable(ReactionVariable, BinaryVariable):
    """
    Class to represent a type of binary variable used to tell whether the
    enzyme is being regulated or not such that:
        EUU + EDU + BigM*EUDUSE <= BigM
    i.e. if the enzyme is being up/down regulated, then EUDUSE is 1
    """

    def __init__(self, reaction, **kwargs):
        if not 'lb' in kwargs:
            kwargs['lb'] = 0
        if not 'ub' in kwargs:
            kwargs['ub'] = 1

        ReactionVariable.__init__(self, reaction,
                                  type=get_binary_type(),
                                  **kwargs)

    prefix = 'EUDUSE_'
