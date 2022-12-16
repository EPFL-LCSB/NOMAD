from pytfa.optim.variables import *

'''
Define Variables for inverseNRA
'''

class TargetDistance(ModelVariable):
    """
    A variable to describe the L2 distance to a target solution
    """
    def __init__(self, model, id_, **kwargs):
        ModelVariable.__init__(self, model, id_, **kwargs)

    prefix = 'TD_'
