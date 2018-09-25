from metapoppy.event import Event
from ..tbcompartments import *


class Replication(Event):
    def __init__(self, cell_type):
        self._cell_type = cell_type
        Event.__init__(self)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._cell_type)

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._cell_type: 1})


class IntracellularBacterialReplication(Replication):
    def __init__(self):
        self._sigmoid = 0
        self._carrying_capacity = 0
        Replication.__init__(self, BACTERIUM_INTRACELLULAR_MACROPHAGE)

    def set_parameters(self, sigmoid, carrying_capacity):
        self._sigmoid = sigmoid
        self._carrying_capacity = carrying_capacity

    def _calculate_state_variable_at_patch(self, network, patch_id):
        bac = network.get_compartment_value(patch_id, BACTERIUM_INTRACELLULAR_MACROPHAGE)
        # Return 0 if no bacteria exist
        if not bac:
            return 0
        mac = network.get_compartment_value(patch_id, MACROPHAGE_INFECTED)
        return bac * (1 - (bac ** self._sigmoid * 1.0 /
                      (bac ** self._sigmoid + (self._carrying_capacity * mac) ** self._sigmoid)))
