from metapoppy import *


class Change(Event):
    def __init__(self, rp_key, compartment_from, compartment_to):
        self._comp_from = compartment_from
        self._comp_to = compartment_to
        Event.__init__(self, rp_key)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._comp_from)

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._comp_from: -1, self._comp_to: 1})


class Infect(Change):
    def __init__(self, rp_key, susceptible_compartment, infectious_compartment, infected_compartment):
        self._infectious = infectious_compartment
        Change.__init__(self, rp_key, susceptible_compartment, infected_compartment)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._comp_from) * \
               network.get_compartment_value(patch_id, self._infectious)
