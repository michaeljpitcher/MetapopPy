from .network import *


class Event(object):

    def __init__(self, dependent_compartments, dependent_attributes, reaction_parameter_key, additional_parameter_keys=None):
        """
        Create an event.
        :param reaction_parameter_key: The parameter key which corresponds to the reaction parameter
        :param additional_parameter_keys: Other parameter keys which are required by this event
        """
        self.dependent_compartments = dependent_compartments
        self.dependent_attributes = dependent_attributes
        self._reaction_parameter_key = reaction_parameter_key
        self._reaction_parameter = 0.0
        self._parameters = {}
        if additional_parameter_keys:
            for p in additional_parameter_keys:
                self._parameters[p] = 0.0

    def set_parameters(self, parameter_values):
        """
        Given a set of parameter values upon experiment configure, assign these to the event
        :param parameter_values:
        :return:
        """
        assert parameter_values[self._reaction_parameter_key] >= 0.0, "Reaction parameter cannot be negative"
        self._reaction_parameter = parameter_values[self._reaction_parameter_key]
        for p in self._parameters:
            self._parameters[p] = parameter_values[p]

    def calculate_rate_at_patch(self, network, patch_id):
        """
        Calculate the rate of this event at this patch. Calculation is reaction parameter (single value, defined at
        experiment initialisation) multiplied by state variable (calculated during simulation)
        :param network:
        :param patch_id:
        :return:
        """
        return self._reaction_parameter * self._calculate_state_variable_at_patch(network, patch_id)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        """
        Determine the state variable at patch (i.e. number of possible occurrences). Must be overriden as specific to
        each event type.
        :param network:
        :param patch_id:
        :return:
        """
        raise NotImplementedError

    def perform(self, network, patch_id):
        """
        Event is performed at a patch, updating it (and other patches). Must be overriden as specific to each event
        type.
        :param network:
        :param patch_id:
        :return:
        """
        raise NotImplementedError


class PatchTypeEvent(Event):
    """
    An event which can only occur at a specific type of patch. Rate will always be zero if calculated at a patch which
    does not match the patch type.
    """
    def __init__(self, patch_type, dependent_compartments, dependent_attributes, reaction_parameter_key,
                 additional_parameter_keys=None):
        """
        Create a patch type event
        :param patch_type:
        :param reaction_parameter_key:
        :param additional_parameter_keys:
        """
        self._patch_type = patch_type
        Event.__init__(self, dependent_compartments, dependent_attributes, reaction_parameter_key,
                       additional_parameter_keys)

    def calculate_rate_at_patch(self, network, patch_id):
        """
        Calculate rate. Zero if at the wrong patch type, otherwise, same as Event.
        :param network:
        :param patch_id:
        :return:
        """
        if network.node[patch_id][TypedNetwork.PATCH_TYPE] == self._patch_type:
            return Event.calculate_rate_at_patch(self, network, patch_id)
        else:
            return 0.0

    def _calculate_state_variable_at_patch(self, patch_data, edges):
        """
        Determine the state variable at patch (i.e. number of possible occurrences). Must be overriden as specific to
        each event type.
        :param network:
        :param patch_id:
        :return:
        """
        raise NotImplementedError

    def perform(self, network, patch_id):
        """
        Event is performed at a patch, updating it (and other patches). Must be overriden as specific to each event
        type.
        :param network:
        :param patch_id:
        :return:
        """
        raise NotImplementedError
