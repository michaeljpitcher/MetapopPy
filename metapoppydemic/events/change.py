from metapoppy import *


class Change(Event):

    RATE_OF_CHANGE = 'rate_of_change_'

    def __init__(self, compartment_from, compartment_to):
        self._comp_from = compartment_from
        self._comp_to = compartment_to
        Event.__init__(self, [compartment_from], [])

    def _define_parameter_keys(self):
        return Change.RATE_OF_CHANGE + self._comp_from + '_' + self._comp_to, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._comp_from)

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._comp_from: -1, self._comp_to: 1})


class Infect(Change):
    INFECTION_RATE_KEY = 'infection_rate'

    def __init__(self, susceptible_compartment, infectious_compartment, infected_compartment):
        self._infectious = infectious_compartment
        Change.__init__(self, susceptible_compartment, infected_compartment)

    def _define_parameter_keys(self):
        return Infect.INFECTION_RATE_KEY, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._comp_from) * \
               network.get_compartment_value(patch_id, self._infectious)
