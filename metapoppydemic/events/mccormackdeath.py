from metapoppy import Event
from ..environment.mccormackenvironment import *


class McCormackDeath(Event):

    RATE_OF_DEATH = 'rate_of_death_'

    def __init__(self, compartment, all_comps):
        self._comp = compartment
        self._all_comps = all_comps
        Event.__init__(self, all_comps, [], [])

    def _define_parameter_keys(self):
        return McCormackDeath.RATE_OF_DEATH + self._comp, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        d1 = network.get_attribute_value(patch_id, McCormackEnvironment.BASE_DEATH_RATE)
        d2 = network.get_attribute_value(patch_id, McCormackEnvironment.POPULATION_DEATH_RATE)
        d = (d1 + (d2 * network.get_compartment_value(patch_id, self._all_comps)))
        return network.get_compartment_value(patch_id, self._comp)* d

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._comp: -1})


class McCormackDeathInfection(Event):

    RATE_OF_DEATH = 'rate_of_infection_death_'

    def __init__(self, compartment):
        self._comp = compartment
        Event.__init__(self, [compartment], [], [])

    def _define_parameter_keys(self):
        return McCormackDeathInfection.RATE_OF_DEATH + self._comp, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        alpha = network.get_attribute_value(patch_id, McCormackEnvironment.INFECTION_DEATH_RATE)
        return alpha * network.get_compartment_value(patch_id, self._comp)

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._comp: -1})
