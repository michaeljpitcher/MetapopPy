from metapoppy import *


class Create(Event):

    RATE_OF_BIRTH = 'rate_of_birth_'

    def __init__(self, compartment):
        self._comp = compartment
        Event.__init__(self, [], [], [])

    def _define_parameter_keys(self):
        return Create.RATE_OF_BIRTH + self._comp, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return 1

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._comp: 1})
