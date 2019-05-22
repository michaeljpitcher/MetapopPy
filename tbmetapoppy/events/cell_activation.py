from ..tbpulmonaryenvironment import TBPulmonaryEnvironment
from metapoppy.event import Event
from parameters import RATE, HALF_SAT


class CellActivation(Event):
    ACTIVATION_BY = '_activation_by_'

    def __init__(self, resting_cell, triggers):
        self._half_sat_key = None
        self._resting_cell = resting_cell
        self._activated_cell = TBPulmonaryEnvironment.ACTIVATED_CELL[self._resting_cell]
        self._triggers = triggers
        Event.__init__(self, [self._resting_cell] + self._triggers, [], [])

    def _define_parameter_keys(self):
        rp_key = self._resting_cell + CellActivation.ACTIVATION_BY + '_'.join(self._triggers) + RATE
        self._half_sat_key = self._resting_cell + CellActivation.ACTIVATION_BY + '_'.join(self._triggers) + HALF_SAT
        return rp_key, [self._half_sat_key]

    def _calculate_state_variable_at_patch(self, network, patch_id):
        trigger_count = network.get_compartment_value(patch_id, self._triggers)
        if not trigger_count:
            return 0
        return network.get_compartment_value(patch_id, self._resting_cell) * \
            (float(trigger_count) / (trigger_count + self._parameters[self._half_sat_key]))

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._resting_cell: -1, self._activated_cell: 1})
