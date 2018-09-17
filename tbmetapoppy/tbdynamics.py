from metapoppy.dynamics import Dynamics
from pulmonarynetwork import *


class TBDynamics(Dynamics):

    BACTERIUM_EXTRACELLULAR_REPLICATING = 'b_er'
    BACTERIUM_EXTRACELLULAR_DORMANT = 'b_ed'
    EXTRACELLULAR_BACTERIA = [BACTERIUM_EXTRACELLULAR_REPLICATING + BACTERIUM_EXTRACELLULAR_DORMANT]
    INTRACELLULAR_BACTERIA = []
    BACTERIA = EXTRACELLULAR_BACTERIA + INTRACELLULAR_BACTERIA

    MACROPHAGES = []

    DENDRITIC_CELLS = []

    T_CELLS = []

    COMPARTMENTS = BACTERIA + MACROPHAGES + DENDRITIC_CELLS + T_CELLS

    def __init__(self):
        # Build network
        pulmonary_network = PulmonaryNetwork(TBDynamics.COMPARTMENTS)
        pulmonary_network.build()
        Dynamics.__init__(self, pulmonary_network)

    def _create_events(self):
        # TODO
        pass

    def _seed_events(self, params):
        # TODO
        pass

    def _seed_network(self, params):
        # TODO
        pass

    def _get_results(self):
        # TODO
        pass
