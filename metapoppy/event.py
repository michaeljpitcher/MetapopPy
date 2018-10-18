from .network import *


class Event(object):

    def __init__(self, reaction_parameter_key, additional_parameter_keys=None):
        self._reaction_parameter_key = reaction_parameter_key
        self._reaction_parameter = 0.0
        self._parameters = {}
        if additional_parameter_keys:
            for p in additional_parameter_keys:
                self._parameters[p] = 0.0

    def set_parameters(self, parameter_values):
        assert parameter_values[self._reaction_parameter_key] >= 0.0, "Reaction parameter cannot be negative"
        self._reaction_parameter = parameter_values[self._reaction_parameter_key]
        for p in self._parameters:
            self._parameters[p] = parameter_values[p]

    def calculate_rate_at_patch(self, network, patch_id):
        return self._reaction_parameter * self._calculate_state_variable_at_patch(network, patch_id)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        raise NotImplementedError

    def perform(self, network, patch_id):
        raise NotImplementedError


class PatchTypeEvent(Event):
    def __init__(self, patch_type, reaction_parameter_key, additional_parameter_keys=None):
        self._patch_type = patch_type
        Event.__init__(self, reaction_parameter_key, additional_parameter_keys)

    def calculate_rate_at_patch(self, network, patch_id):
        if network.node[patch_id][TypedNetwork.PATCH_TYPE] == self._patch_type:
            return Event.calculate_rate_at_patch(self, network, patch_id)
        else:
            return 0.0

    def _calculate_state_variable_at_patch(self, patch_data, edges):
        raise NotImplementedError

    def perform(self, network, patch_id):
        raise NotImplementedError
