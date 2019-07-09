from metapoppy import Event
from ..environment.mccormackenvironment import *


class McCormackInfection(Event):

    RATE_OF_INFECTION = 'rate_of_infection_'

    def __init__(self, compartment_s, compartment_i, all_comps):
        self._comp_s = compartment_s
        self._comp_i = compartment_i
        self._all_comps = all_comps
        Event.__init__(self, all_comps, [], [])

    def _define_parameter_keys(self):
        return McCormackInfection.RATE_OF_INFECTION + self._comp_s + '_' + self._comp_i, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        lambda_i = network.get_attribute_value(patch_id, McCormackEnvironment.INFECTION_LAMBDA)
        return (lambda_i / network.get_compartment_value(patch_id, self._all_comps)) * \
               network.get_compartment_value(patch_id, self._comp_s) * \
               network.get_compartment_value(patch_id, self._comp_i)

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._comp_s: -1, self._comp_i:1})
