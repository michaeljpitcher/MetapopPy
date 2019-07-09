from metapoppy import Event
from ..environment.mccormackenvironment import *


class McCormackRecover(Event):

    RATE_OF_RECOVERY = 'rate_of_recovery_'

    def __init__(self, compartment_from, compartment_to):
        self._comp_from = compartment_from
        self._comp_to = compartment_to
        Event.__init__(self, [compartment_from], [], [])

    def _define_parameter_keys(self):
        return McCormackRecover.RATE_OF_RECOVERY + self._comp_from, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        alpha = network.get_attribute_value(patch_id, McCormackEnvironment.RECOVERY_RATE)
        return alpha * network.get_compartment_value(patch_id, self._comp_from)

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._comp_from: -1, self._comp_to: 1})
