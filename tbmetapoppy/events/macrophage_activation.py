from ..tbcompartments import *
from metapoppy.event import Event


class MacrophageActivationByTCell(Event):

    HALF_SAT = 0

    def __init__(self):
        Event.__init__(self)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        t_cell = network.get_compartment_value(patch_id, T_CELL_ACTIVATED)
        if not t_cell:
            return 0
        return network.get_compartment_value(patch_id, MACROPHAGE_RESTING) * \
            (t_cell / (t_cell + MacrophageActivationByTCell.HALF_SAT))

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {MACROPHAGE_RESTING: -1, MACROPHAGE_ACTIVATED: 1})


class MacrophageActivationByBacteria(Event):
    HALF_SAT = 0

    def __init__(self):
        Event.__init__(self)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        bac = network.get_compartment_value(patch_id, BACTERIUM_EXTRACELLULAR_REPLICATING) + \
              network.get_compartment_value(patch_id, BACTERIUM_EXTRACELLULAR_DORMANT)
        if not bac:
            return 0
        return network.get_compartment_value(patch_id, MACROPHAGE_RESTING) * \
               (bac / (bac + MacrophageActivationByBacteria.HALF_SAT))

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {MACROPHAGE_RESTING: -1, MACROPHAGE_ACTIVATED: 1})

