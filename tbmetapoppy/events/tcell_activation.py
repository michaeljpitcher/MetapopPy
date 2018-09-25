from tbmetapoppy.pulmonarynetwork import PulmonaryNetwork
from ..tbcompartments import *
from metapoppy.event import PatchTypeEvent


class TCellActivation(PatchTypeEvent):
    def __init__(self, antigen_presenter):
        self._antigen_presenter = antigen_presenter
        self._half_sat = 0
        PatchTypeEvent.__init__(self, PulmonaryNetwork.LYMPH_PATCH)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        apc_count = network.get_compartment_value(patch_id, self._antigen_presenter)
        if not apc_count:
            return 0
        return network.get_compartment_value(patch_id, T_CELL_NAIVE) * \
               (apc_count / (apc_count + self._half_sat))

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {T_CELL_NAIVE: -1, T_CELL_ACTIVATED: 1})
