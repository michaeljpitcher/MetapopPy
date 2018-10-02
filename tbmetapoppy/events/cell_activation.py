from ..tbcompartments import *
from metapoppy.event import Event


class CellActivation(Event):
    def __init__(self, resting_cell, triggers):
        self._half_sat = 0
        self._resting_cell = resting_cell
        self._activated_cell = ACTIVATED_CELL[self._resting_cell]
        self._triggers = triggers
        Event.__init__(self)

    def set_parameters(self, reaction_parameter, half_sat):
        self.set_reaction_parameter(reaction_parameter)
        self._half_sat = half_sat

    def _calculate_state_variable_at_patch(self, network, patch_id):
        trigger_count = sum([network.get_compartment_value(patch_id, n) for n in self._triggers])
        if not trigger_count:
            return 0
        return network.get_compartment_value(patch_id, self._resting_cell) * \
            (float(trigger_count) / (trigger_count + self._half_sat))

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._resting_cell: -1, self._activated_cell: 1})
