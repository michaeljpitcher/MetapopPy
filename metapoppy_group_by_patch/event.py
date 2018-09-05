class Event(object):
    """
    A MetapopPy event. A class of event has a reaction parameter. Each instance of the event class occurs at a single
    patch, and calculates a state variable based on the composition of the subpopulation at that patch, the
    environmental attributes of the patch and the network edges originating from the patch. When performed the event
    updates the patch (and possibly it's neighbours) in some manner.
    """

    # Reaction parameter should be set at sub-class level.
    REACTION_PARAMETER = 0.0

    def __init__(self):
        """
        Create an event. Default is unattached to any patch.
        """
        self._rate = 0
        self._patch_id = None

    def set_patch_id(self, patch_id):
        """
        Event has been assigned a patch, so update the patch ID.
        :param patch_id:
        :return:
        """
        self._patch_id = patch_id

    def rate(self):
        """
        Get function for rate
        :return:
        """
        return self._rate

    def update_rate(self, patch, edges):
        """
        The patch and/or edges have been altered - update the event's state variable
        :param patch:
        :param edges:
        :return:
        """
        sv = self._calculate_state_variable(patch, edges)
        self._rate = self.__class__.REACTION_PARAMETER * sv

    def _calculate_state_variable(self, patch, edges):
        """
        Calculate the state variable that determines rate of this event occuring. Must be overridden by sub-class.
        :param patch:
        :param edges:
        :return:
        """
        raise NotImplementedError

    def perform(self, network):
        """
        Event is performed and updates the network. Must be overridden by sub-class.
        :param network:
        :return:
        """
        raise NotImplementedError
