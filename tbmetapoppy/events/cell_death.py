from metapoppy.event import Event
from ..tbcompartments import *


class CellDeath(Event):
    def __init__(self, dying_compartment):
        self._dying_compartment = dying_compartment
        if dying_compartment == MACROPHAGE_INFECTED:
            self._internal_compartment = BACTERIUM_INTRACELLULAR_MACROPHAGE
        elif dying_compartment == DENDRITIC_CELL_MATURE:
            self._internal_compartment = BACTERIUM_INTRACELLULAR_DENDRITIC
        else:
            self._internal_compartment = None
        Event.__init__(self)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._dying_compartment)

    def perform(self, network, patch_id):
        changes = {self._dying_compartment: -1}
        if self._internal_compartment:
            bac_to_release = int(round(float(
                network.get_compartment_value(patch_id, self._internal_compartment)) /
                                       network.get_compartment_value(patch_id, self._dying_compartment)))
            changes[self._internal_compartment] = bac_to_release * -1
            changes[BACTERIUM_EXTRACELLULAR_DORMANT] = bac_to_release
        network.update_patch(patch_id, changes)


class MacrophageBursting(CellDeath):
    def __init__(self):
        self._sigmoid = 0
        self._capacity = 0
        CellDeath.__init__(self, MACROPHAGE_INFECTED)

    def set_parameters(self, sigmoid, capacity):
        self._sigmoid = sigmoid
        self._capacity = capacity

    def _calculate_state_variable_at_patch(self, network, patch_id):
        bac = network.get_compartment_value(patch_id, BACTERIUM_INTRACELLULAR_MACROPHAGE)
        if not bac:
            return 0
        mac = network.get_compartment_value(patch_id, self._dying_compartment)
        return mac * ((float(bac) ** self._sigmoid) / (bac ** self._sigmoid +
                ((self._capacity * mac) ** self._sigmoid)))


class TCellDestroysMacrophage(CellDeath):
    def __init__(self):
        self._half_sat = 0
        CellDeath.__init__(self, MACROPHAGE_INFECTED)

    def set_parameters(self, half_sat):
        self._half_sat = half_sat

    def _calculate_state_variable_at_patch(self, network, patch_id):
        tcell = network.get_compartment_value(patch_id, T_CELL_ACTIVATED)
        if not tcell:
            return 0
        mac = network.get_compartment_value(patch_id, self._dying_compartment)
        return mac * (float(tcell) / (tcell + self._half_sat))


class MacrophageDestroysBacterium(CellDeath):
    def __init__(self, macrophage_type, bacterium_type):
        self._macrophage_type = macrophage_type
        self._bacterium_type = bacterium_type
        self._half_sat = 0
        CellDeath.__init__(self, bacterium_type)

    def set_parameters(self, half_sat):
        self._half_sat = half_sat

    def _calculate_state_variable_at_patch(self, network, patch_id):
        mac = network.get_compartment_value(patch_id, self._macrophage_type)
        if not mac:
            return 0
        return network.get_compartment_value(patch_id, self._bacterium_type) * float(mac) / (mac + self._half_sat)
