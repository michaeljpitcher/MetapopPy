from ..tbpulmonaryenvironment import TBPulmonaryEnvironment
from metapoppy.event import PatchTypeEvent
from parameters import RATE, SIGMOID, HALF_SAT


class BacteriumChangeStateThroughOxygen(PatchTypeEvent):

    BACTERIUM_CHANGE_TO_DORMANT = 'bacterium_change_to_dormant'
    BACTERIUM_CHANGE_TO_REPLICATING = 'bacterium_change_to_replicating'

    def __init__(self, replicating_to_dormant):
        if replicating_to_dormant:
            self._compartment_from = TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_REPLICATING
            self._compartment_to = TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT
            self._key_prefix = BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_DORMANT
        else:
            self._compartment_from = TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT
            self._compartment_to = TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_REPLICATING
            self._key_prefix = BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_REPLICATING
        self._sigmoid_key = None
        self._half_sat_key = None
        PatchTypeEvent.__init__(self, TBPulmonaryEnvironment.ALVEOLAR_PATCH, [self._compartment_from],
                                [TBPulmonaryEnvironment.OXYGEN_TENSION], [])

    def _define_parameter_keys(self):
        self._sigmoid_key = self._key_prefix + SIGMOID
        self._half_sat_key = self._key_prefix + HALF_SAT
        return self._key_prefix + RATE, [self._sigmoid_key, self._half_sat_key]

    def _calculate_state_variable_at_patch(self, network, patch_id):
        bac = network.get_compartment_value(patch_id, self._compartment_from)
        if not bac:
            return 0
        sig = self._parameters[self._sigmoid_key]
        halfsat = self._parameters[self._half_sat_key]
        o2 = network.get_attribute_value(patch_id, TBPulmonaryEnvironment.OXYGEN_TENSION)
        if o2 == 0:
            return bac
        return bac * ((o2 ** sig) / (halfsat ** sig + o2 ** sig))

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._compartment_from: -1, self._compartment_to: 1})
