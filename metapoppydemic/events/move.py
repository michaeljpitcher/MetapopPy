from metapoppy import *
import numpy


class Move(Event):
    def __init__(self, rp_key, moving_compartment):
        self._mover = moving_compartment
        Event.__init__(self, rp_key)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._mover) * (len(network.edges([patch_id])) > 0)

    def perform(self, network, patch_id):
        # Moves along an edge at random
        edges = network.edges([patch_id], data=True)
        chosen_neighbour, _ = numpy.random.choice(edges)
        network.update_patch(patch_id, {self._mover: -1})
        network.update_patch(chosen_neighbour, {self._mover: 1})
