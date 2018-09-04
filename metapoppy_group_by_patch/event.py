class Event(object):

    REACTION_PARAMETER = 0.0

    def __init__(self):
        self._rate = 0
        self._patch_id = None

    def set_patch_id(self, patch_id):
        self._patch_id = patch_id

    def rate(self):
        return self._rate

    def update_rate(self, patch, edges):
        sv = self._calculate_state_variable(patch, edges)
        self._rate = self.__class__.REACTION_PARAMETER * sv

    def _calculate_state_variable(self, patch, edges):
        raise NotImplementedError

    def perform(self, network):
        raise NotImplementedError
