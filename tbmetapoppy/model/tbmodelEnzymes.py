from tbmodel import *
from ..events import *


class TBModelEnzymes(TBDynamics):

    def __init__(self, network_config):
        TBDynamics.__init__(self, network_config)

    def _create_events(self):
        events = TBDynamics._create_events(self)

        cas_liq = CaseumLiquefy()

        events.append(cas_liq)
        return events