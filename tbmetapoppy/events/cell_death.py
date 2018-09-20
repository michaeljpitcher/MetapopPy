from tbmetapoppy import *
from metapoppy.event import Event


class CellDeath(Event):
    def __init__(self, cell_type):
        self.cell_type = cell_type
        Event.__init__(self)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self.cell_type)

    def perform(self, network, patch_id):
        changes = {self.cell_type: -1}
        if self.cell_type == TBDynamics.MACROPHAGE_INFECTED:
            bac_to_release = int(round(float(
                network.get_compartment_value(patch_id, TBDynamics.BACTERIUM_INTRACELLULAR_MACROPHAGE)) /
                 network.get_compartment_value(patch_id, TBDynamics.MACROPHAGE_INFECTED)))
            changes[TBDynamics.BACTERIUM_INTRACELLULAR_MACROPHAGE] = bac_to_release * -1
            changes[TBDynamics.BACTERIUM_EXTRACELLULAR_DORMANT] = bac_to_release
        elif self.cell_type == TBDynamics.DENDRITIC_CELL_MATURE:
            changes[TBDynamics.BACTERIUM_INTRACELLULAR_DENDRITIC] = -1
            changes[TBDynamics.BACTERIUM_EXTRACELLULAR_DORMANT] = 1
        network.update_patch(patch_id, changes)


class MacrophageBursting(CellDeath):
    SIGMOID = 0
    CARRYING_CAPACITY = 0

    def __init__(self):
        CellDeath.__init__(self, TBDynamics.MACROPHAGE_INFECTED)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        bac = network.get_compartment_value(patch_id, TBDynamics.BACTERIUM_INTRACELLULAR_MACROPHAGE)
        if not bac:
            return 0
        mac = network.get_compartment_value(patch_id, TBDynamics.MACROPHAGE_INFECTED)
        return mac * ((float(bac) ** MacrophageBursting.SIGMOID) / (bac ** MacrophageBursting.SIGMOID +
                ((MacrophageBursting.CARRYING_CAPACITY * mac) ** MacrophageBursting.SIGMOID)))



class MacrophageDestroysBacterium(CellDeath):
    def __init__(self, macrophage_type, bacterium_type):
        self._macrophage_type = macrophage_type
        self._bacterium_type = bacterium_type
        CellDeath.__init__(self, bacterium_type)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._bacterium_type) * \
               network.get_compartment_value(patch_id, self._macrophage_type)


class TCellDestroysMacrophage(CellDeath):
    HALF_SAT = 0

    def __init__(self):
        CellDeath.__init__(self, TBDynamics.MACROPHAGE_INFECTED)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        tcell = network.get_compartment_value(patch_id, TBDynamics.T_CELL_ACTIVATED)
        if not tcell:
            return 0
        mac = network.get_compartment_value(patch_id, TBDynamics.MACROPHAGE_INFECTED)
        return mac * (float(tcell) / (tcell + TCellDestroysMacrophage.HALF_SAT))
