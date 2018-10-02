from tbmetapoppy.pulmonarynetwork import PulmonaryNetwork
from metapoppy.event import PatchTypeEvent
from ..tbcompartments import *

# TODO - more advanced dormancy dynamics


class BacteriumBecomesDormant(PatchTypeEvent):

    def __init__(self):
        self._sigmoid = 0
        self._half_sat = 0
        PatchTypeEvent.__init__(self,  PulmonaryNetwork.ALVEOLAR_PATCH)

    def set_parameters(self, reaction_parameter, sigmoid, half_sat):
        self.set_reaction_parameter(reaction_parameter)
        assert sigmoid < 0, "Sigmoid must be negative"
        self._sigmoid = sigmoid
        self._half_sat = half_sat

    def _calculate_state_variable_at_patch(self, network, patch_id):
        bac = network.get_compartment_value(patch_id, BACTERIUM_EXTRACELLULAR_REPLICATING)
        if not bac:
            return 0
        o2 = network.get_attribute_value(patch_id, PulmonaryNetwork.OXYGEN_TENSION)
        return bac * ((o2 ** self._sigmoid) / (self._half_sat ** self._sigmoid + o2 ** self._sigmoid))

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {BACTERIUM_EXTRACELLULAR_REPLICATING: -1, BACTERIUM_EXTRACELLULAR_DORMANT: 1})


class BacteriumBecomesReplicating(PatchTypeEvent):

    def __init__(self):
        self._sigmoid = 0
        self._half_sat = 0
        PatchTypeEvent.__init__(self, PulmonaryNetwork.ALVEOLAR_PATCH)

    def set_parameters(self, reaction_parameter, sigmoid, half_sat):
        self.set_reaction_parameter(reaction_parameter)
        assert sigmoid > 0, "Sigmoid must be positive"
        self._sigmoid = sigmoid
        self._half_sat = half_sat

    def _calculate_state_variable_at_patch(self, network, patch_id):
        bac = network.get_compartment_value(patch_id, BACTERIUM_EXTRACELLULAR_DORMANT)
        o2 = network.get_attribute_value(patch_id, PulmonaryNetwork.OXYGEN_TENSION)
        return bac * ((o2 ** self._sigmoid) / (self._half_sat ** self._sigmoid + o2 ** self._sigmoid))

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {BACTERIUM_EXTRACELLULAR_DORMANT: -1, BACTERIUM_EXTRACELLULAR_REPLICATING: 1})
