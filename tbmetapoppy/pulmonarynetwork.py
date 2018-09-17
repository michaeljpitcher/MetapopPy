from metapoppy import *


class PulmonaryNetwork(TypedNetwork):

    VENTILATION = 'ventilation'
    PERFUSION = 'perfusion'
    OXYGEN_TENSION = 'oxygen_tension'
    DRAINAGE = 'drainage'
    PATCH_ATTRIBUTES = [VENTILATION, PERFUSION, OXYGEN_TENSION, DRAINAGE]

    # TODO - more edge attributes
    EDGE_ATTRIBUTES = [PERFUSION]

    def __init__(self, compartments):
        TypedNetwork.__init__(self, compartments, PulmonaryNetwork.PATCH_ATTRIBUTES, PulmonaryNetwork.EDGE_ATTRIBUTES)

    def build(self):
        # TODO - build a space-filling tree
        pass
