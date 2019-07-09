from metapoppy import *
import numpy
from ..environment.mccormackenvironment import *


class McCormackMove(Event):
    MOVEMENT_RATE_KEY = 'movement_rate_'

    def __init__(self, moving_compartment):
        self._mover = moving_compartment
        Event.__init__(self, [moving_compartment], [], [])

    def _define_parameter_keys(self):
        return McCormackMove.MOVEMENT_RATE_KEY + self._mover, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        k = network.get_attribute_value(patch_id, McCormackEnvironment.CARRYING_CAPACITY)
        return (1.0/k) * network.get_compartment_value(patch_id, self._mover) * len(network.edges([patch_id]))

    def perform(self, network, patch_id):
        # Moves along an edge at random
        edges = [v for _,v,_ in network.edges([patch_id], data=True)]
        chosen_neighbour = numpy.random.choice(edges)
        network.update_patch(patch_id, {self._mover: -1})
        network.update_patch(chosen_neighbour, {self._mover: 1})
