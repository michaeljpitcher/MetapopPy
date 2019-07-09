from metapoppy import *


class Death(Event):

    RATE_OF_DEATH = 'rate_of_death_'

    def __init__(self, compartment):
        self._comp = compartment
        Event.__init__(self, [compartment], [], [])

    def _define_parameter_keys(self):
        return Death.RATE_OF_DEATH + self._comp, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._comp)

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._comp: -1})
