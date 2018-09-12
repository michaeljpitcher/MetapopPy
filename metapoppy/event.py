from .network import *


class Event(object):

    def __init__(self):
        self._reaction_parameter = 0.0
        self._handler = None

    def set_reaction_parameter(self, value):
        self._reaction_parameter = value

    def calculate_rate_at_patch(self, network, patch_id):
        return self._reaction_parameter * self._calculate_state_variable_at_patch(network, patch_id)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        raise NotImplementedError

    def perform(self, network, patch_id):
        raise NotImplementedError


class PatchTypeEvent(Event):
    def __init__(self, patch_type):
        self._patch_type = patch_type
        Event.__init__(self)

    def calculate_rate_at_patch(self, network, patch_id):
        patch_data = network.node[patch_id]
        if patch_data[TypedNetwork.PATCH_TYPE] == self._patch_type:
            return self._reaction_parameter * self._calculate_state_variable_at_patch(network, patch_id)
        else:
            return 0.0

    def _calculate_state_variable_at_patch(self, patch_data, edges):
        raise NotImplementedError

    def perform(self, network, patch_id):
        raise NotImplementedError
