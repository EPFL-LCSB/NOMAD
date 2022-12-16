from pytfa.optim.constraints import *
'''
Define new constraints for inverseNRA
'''

class TargetDistanceCoupling(GenericConstraint):
    """
    Class to couple the L2 distance for a target solution to a distance variable
    D - SUM((logv_i - logv_target_i)**2) = 0
    """
    def __init__(self, model, expr, id_, **kwargs):

        GenericConstraint.__init__(self,
                                   id_=id_,
                                   expr=expr,
                                   model=model,
                                   **kwargs)

    prefix = 'TDC_'
