from tbmetapoppy import *
from metapoppy.event import PatchTypeEvent


class CellRecruitmentLung(PatchTypeEvent):
    def __init__(self, cell_type):
        self._cell_type = cell_type
        PatchTypeEvent.__init__(self, PulmonaryNetwork.ALVEOLAR_PATCH)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        # TODO - build in enhanced recruitment levels
        return network.get_attribute_value(patch_id, PulmonaryNetwork.PERFUSION)

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._cell_type: 1})


class CellRecruitmentLymph(PatchTypeEvent):
    def __init__(self, cell_type):
        self._cell_type = cell_type
        PatchTypeEvent.__init__(self, PulmonaryNetwork.LYMPH_PATCH)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        # TODO - build in enhanced recruitment levels
        return 1

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._cell_type: 1})
