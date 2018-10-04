from metapoppy.event import Event
from ..tbcompartments import *


class Replication(Event):
    def __init__(self, replication_rate_key, cell_type, additional_parameter_keys=None):
        self._cell_type = cell_type
        Event.__init__(self, replication_rate_key, additional_parameter_keys)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._cell_type)

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._cell_type: 1})


class IntracellularBacterialReplication(Replication):
    def __init__(self, replication_rate_key, sigmoid_key, capacity_key):
        self._sigmoid_key = sigmoid_key
        self._carrying_capacity_key = capacity_key
        Replication.__init__(self, replication_rate_key, BACTERIUM_INTRACELLULAR_MACROPHAGE,
                             [sigmoid_key, capacity_key])

    def _calculate_state_variable_at_patch(self, network, patch_id):
        bac = network.get_compartment_value(patch_id, BACTERIUM_INTRACELLULAR_MACROPHAGE)
        # Return 0 if no bacteria exist
        if not bac:
            return 0
        sig = self._parameters[self._sigmoid_key]
        cap = self._parameters[self._carrying_capacity_key]
        mac = network.get_compartment_value(patch_id, MACROPHAGE_INFECTED)
        return bac * (1 - (float(bac ** sig) / (bac ** sig + (cap * mac) ** sig)))
