from ..tbpulmonarynetwork import TBPulmonaryNetwork
from metapoppy.event import PatchTypeEvent
# TODO - more advanced dormancy dynamics


class BacteriumChangeStateThroughOxygen(PatchTypeEvent):

    def __init__(self, replicating_to_dormant, change_rate_key, sigmoid_key, half_sat_key):
        if replicating_to_dormant:
            self._compartment_from = TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING
            self._compartment_to = TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT
        else:
            self._compartment_from = TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT
            self._compartment_to = TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING
        self._sigmoid_key = sigmoid_key
        self._half_sat_key = half_sat_key
        PatchTypeEvent.__init__(self, TBPulmonaryNetwork.ALVEOLAR_PATCH, [self._compartment_from],
                                [TBPulmonaryNetwork.OXYGEN_TENSION], change_rate_key,
                                [sigmoid_key, half_sat_key])

    def _calculate_state_variable_at_patch(self, network, patch_id):
        bac = network.get_compartment_value(patch_id, self._compartment_from)
        if not bac:
            return 0
        sig = self._parameters[self._sigmoid_key]
        halfsat = self._parameters[self._half_sat_key]
        o2 = network.get_attribute_value(patch_id, TBPulmonaryNetwork.OXYGEN_TENSION)
        return bac * ((o2 ** sig) / (halfsat ** sig + o2 ** sig))

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._compartment_from: -1, self._compartment_to: 1})
