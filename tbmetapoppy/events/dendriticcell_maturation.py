from tbmetapoppy import *
from metapoppy.event import Event


class DendriticCellMaturation(Event):
    def __init__(self, bacterium_type):
        self._bacterium_type = bacterium_type
        Event.__init__(self)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, TBDynamics.DENDRITIC_CELL_IMMATURE) * \
               network.get_compartment_value(patch_id, self._bacterium_type)

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {TBDynamics.DENDRITIC_CELL_IMMATURE: -1,
                                        TBDynamics.DENDRITIC_CELL_MATURE: 1,
                                        self._bacterium_type: -1,
                                        TBDynamics.BACTERIUM_INTRACELLULAR_DENDRITIC: 1})
