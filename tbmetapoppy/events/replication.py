from metapoppy.event import *
from ..tbpulmonaryenvironment import *
from parameters import *


class Replication(Event):

    REPLICATION_RATE = '_replication_rate'

    def __init__(self, cell_type):
        self._cell_type = cell_type
        Event.__init__(self, [self._cell_type], [], [])

    def _define_parameter_keys(self):
        return self._cell_type + Replication.REPLICATION_RATE, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._cell_type)

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._cell_type: 1})


class IntracellularBacterialReplication(Replication):

    def __init__(self):
        Replication.__init__(self, TBPulmonaryEnvironment.BACTERIUM_INTRACELLULAR_MACROPHAGE)
        self._dependent_compartments += [TBPulmonaryEnvironment.BACTERIUM_INTRACELLULAR_MACROPHAGE,
                                         TBPulmonaryEnvironment.MACROPHAGE_INFECTED]

    def _define_parameter_keys(self):
        return self._cell_type + Replication.REPLICATION_RATE, [INTRACELLULAR_REPLICATION_SIGMOID,
                                                                MACROPHAGE_CAPACITY]

    def _calculate_state_variable_at_patch(self, network, patch_id):
        bac = network.get_compartment_value(patch_id, TBPulmonaryEnvironment.BACTERIUM_INTRACELLULAR_MACROPHAGE)
        # Return 0 if no bacteria exist
        if not bac:
            return 0
        sig = self._parameters[INTRACELLULAR_REPLICATION_SIGMOID]
        cap = self._parameters[MACROPHAGE_CAPACITY]
        mac = network.get_compartment_value(patch_id, TBPulmonaryEnvironment.MACROPHAGE_INFECTED)
        return bac * (1 - (float(bac ** sig) / (bac ** sig + (cap * mac) ** sig)))


class PatchTypeReplication(PatchTypeEvent):

    REPLICATION_RATE = '_replication_rate'
    CAPACITY = '_replication_capacity'
    SIGMOID = '_replication_sigmoid'

    def __init__(self, patch_type, cell_type):
        self._cell_type = cell_type
        PatchTypeEvent.__init__(self, patch_type, [self._cell_type], [], [])

    def _define_parameter_keys(self):
        self.cap_key = self._cell_type + '_' + self._patch_type + PatchTypeReplication.CAPACITY
        self.sig_key = self._cell_type + '_' + self._patch_type + PatchTypeReplication.SIGMOID
        return self._cell_type + '_' + self._patch_type + PatchTypeReplication.REPLICATION_RATE, [self.cap_key, self.sig_key]

    def _calculate_state_variable_at_patch(self, network, patch_id):
        cell = network.get_compartment_value(patch_id, self._cell_type)
        if not cell:
            return 0
        sig = self._parameters[self.sig_key]
        return cell * (1 - float(cell)**sig / (cell**sig + self._parameters[self.cap_key]**sig))

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._cell_type: 1})
