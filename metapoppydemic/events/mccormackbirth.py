from metapoppy import Event
from ..environment.mccormackenvironment import *


class McCormackBirth(Event):

    RATE_OF_BIRTH = 'rate_of_birth_'

    def __init__(self, compartment, all_comps):
        self._comp = compartment
        self._all_comps = all_comps
        Event.__init__(self, all_comps, [], [])

    def _define_parameter_keys(self):
        return McCormackBirth.RATE_OF_BIRTH + self._comp, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        b = network.get_attribute_value(patch_id, McCormackEnvironment.BIRTH_RATE)
        return b * network.get_compartment_value(patch_id, self._all_comps)

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._comp: 1})
