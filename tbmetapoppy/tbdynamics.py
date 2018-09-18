from metapoppy.dynamics import Dynamics
from pulmonarynetwork import *
from events import *


class TBDynamics(Dynamics):

    BACTERIUM_EXTRACELLULAR_REPLICATING = 'b_er'
    BACTERIUM_EXTRACELLULAR_DORMANT = 'b_ed'
    EXTRACELLULAR_BACTERIA = [BACTERIUM_EXTRACELLULAR_REPLICATING, BACTERIUM_EXTRACELLULAR_DORMANT]
    INTRACELLULAR_BACTERIA = []
    BACTERIA = EXTRACELLULAR_BACTERIA + INTRACELLULAR_BACTERIA

    MACROPHAGES = []

    DENDRITIC_CELLS = []

    T_CELLS = []

    TB_COMPARTMENTS = BACTERIA + MACROPHAGES + DENDRITIC_CELLS + T_CELLS

    def __init__(self, network_config):
        # Build network
        pulmonary_network = PulmonaryNetwork(network_config, TBDynamics.TB_COMPARTMENTS)
        Dynamics.__init__(self, pulmonary_network)

    def _create_events(self):
        events = []
        for b in TBDynamics.EXTRACELLULAR_BACTERIA:
            events.append(ExtracellularBacterialReplication(b))
        return events

    def _seed_events(self, params):
        # TODO
        pass

    def _seed_network(self, params):
        # TODO
        pass

    def _get_results(self):
        # TODO
        pass
