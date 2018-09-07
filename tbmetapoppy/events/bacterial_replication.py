from metapoppy_group_by_patch.event import Event
from ..pulmonary_network import PulmonaryNetwork


class ExtracellularBacterialReplication(Event):
    def __init__(self, bacterial_type):
        assert bacterial_type in PulmonaryNetwork.EXTRACELLULAR_BACTERIA, "Bacterium type must be extracellular (" \
                                                                          "dormant or replicating)"
        self._bacterial_type = bacterial_type
        Event.__init__(self)

    def _calculate_state_variable(self, patch, edges):
        return patch[self._bacterial_type]

    def perform(self, network):
        network.update_patch(self._patch_id, {self._bacterial_type: 1})
