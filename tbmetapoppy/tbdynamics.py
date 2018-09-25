from metapoppy.dynamics import Dynamics
from .pulmonarynetwork import *
from tbmetapoppy.events import *
from tbcompartments import *


class TBDynamics(Dynamics):

    # Event Parameter keys
    RP_REPLICATION_BER = 'rate_of_extracellular_replication_replicating'
    RP_REPLICATION_BED = 'rate_of_extracellular_replication_dormant'
    RP_MAC_RECRUIT_LUNG = 'rate_of_macrophage_recruitment_lung'

    def __init__(self, network_config):
        # Build network
        pulmonary_network = PulmonaryNetwork(network_config, TB_COMPARTMENTS)
        Dynamics.__init__(self, pulmonary_network)

    def _create_events(self):
        events = []
        # Bacteria replication
        self._ber_replication = Replication(BACTERIUM_EXTRACELLULAR_REPLICATING)
        events.append(self._ber_replication)
        self._bed_replication = Replication(BACTERIUM_EXTRACELLULAR_DORMANT)
        events.append(self._bed_replication)

        # Macrophage recruitment

        return events

    def _seed_events(self, params):
        self._ber_replication.set_reaction_parameter(params[TBDynamics.RP_REPLICATION_BER])
        self._bed_replication.set_reaction_parameter(params[TBDynamics.RP_REPLICATION_BED])

    def _seed_network(self, params):
        # Seed attributes
        self._network.seed_pulmonary_attributes()
        # Seed initial compartments
        # TODO

    def _get_results(self):
        # TODO
        pass
