from tbmetapoppy import *
from metapoppy.event import Event


class MacrophageInfection(Event):

    HALF_SAT = 0

    def __init__(self, bac_type):
        self._bacteria = bac_type
        Event.__init__(self)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        counts = {TBDynamics.BACTERIUM_EXTRACELLULAR_REPLICATING:
                      network.get_compartment_value(patch_id, TBDynamics.BACTERIUM_EXTRACELLULAR_REPLICATING),
                  TBDynamics.BACTERIUM_EXTRACELLULAR_DORMANT:
                      network.get_compartment_value(patch_id, TBDynamics.BACTERIUM_EXTRACELLULAR_DORMANT)}
        total_bac = sum(counts.values())
        if not total_bac:
            return 0
        return network.get_compartment_value(patch_id, TBDynamics.MACROPHAGE_RESTING) * \
            (counts[self._bacteria] / (total_bac + MacrophageInfection.HALF_SAT))

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {TBDynamics.MACROPHAGE_RESTING: -1, TBDynamics.MACROPHAGE_ACTIVATED: 1})
