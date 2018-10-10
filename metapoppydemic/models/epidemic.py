from metapoppy import *
from ..events import *
from compartments import *


class Epidemic(Dynamics):

    def __init__(self, compartments, template_network):
        g = Network(compartments, [], [], template=template_network)
        Dynamics.__init__(self, g)

    def _create_events(self):
        raise NotImplementedError

    def _seed_prototype_network(self, params):
        raise NotImplementedError
