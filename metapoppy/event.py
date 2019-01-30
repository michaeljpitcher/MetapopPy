from .network import *


class Event(object):
    """
    A MetapopPy event which occurs upon a networked metapopulation.

    Each event contains a reaction parameter (the rate of a single occurrence of the event) and a function to calculate
    the state variable (number of possible occurrences) at a patch. These are multiplied together to give the rate of
    occurrence of the event at a patch.

    Each patch contains a function to alter the network in some manner at a patch when the event occurs.

    Patches must define the compartments and attributes their state variable functions are dependent upon (needed to
    propagate patch updates). They must also define the parameter keys that are required for state variable calculation
    (to be updated when the parameters update).
    """

    def __init__(self, dependent_compartments, dependent_patch_attributes, dependent_edge_attributes):
        """
        Create an event.
        """
        self._dependent_compartments = dependent_compartments
        self._dependent_patch_attributes = dependent_patch_attributes
        self._dependent_edge_attributes = dependent_edge_attributes
        self._reaction_parameter_key, self._parameter_keys = self._define_parameter_keys()
        self._parameters = {}
        if self._parameter_keys:
            self._parameters = {p:0.0 for p in self._parameter_keys}
        self._reaction_parameter = 0.0

    def get_dependent_compartments(self):
        return self._dependent_compartments

    def get_dependent_patch_attributes(self):
        return self._dependent_patch_attributes

    def get_dependent_edge_attributes(self):
        return self._dependent_edge_attributes

    def _define_parameter_keys(self):
        raise NotImplementedError

    def reaction_parameter(self):
        return self._reaction_parameter_key

    def parameter_keys(self):
        return [self._reaction_parameter_key] + self._parameter_keys

    def set_parameters(self, parameter_values):
        """
        Given a set of parameter values upon experiment configure, assign these to the event
        :param parameter_values:
        :return:
        """
        assert self._reaction_parameter_key in parameter_values, "Reaction parameter {0} missing".format(
            self._reaction_parameter_key)
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
    def __init__(self, patch_type, dependent_compartments, dependent_attributes, dependent_edge_attributes):
        self._patch_type = patch_type
        Event.__init__(self, dependent_compartments, dependent_attributes, dependent_edge_attributes)

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
