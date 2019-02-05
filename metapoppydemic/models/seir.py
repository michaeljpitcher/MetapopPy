from metapoppy import *
from epidemic import *
from ..events import *
from compartments import *


class SEIRDynamics(Epidemic):

    INIT_S = 'initial_population_susceptible'
    INIT_I = 'initial_population_infected'

    def __init__(self, template_network):
        Epidemic.__init__(self, [SUSCEPTIBLE, EXPOSED, INFECTIOUS, RECOVERED], template_network)

    def _create_events(self):
        # Contact with I moves S to E
        infect = Infect(SUSCEPTIBLE, INFECTIOUS, EXPOSED)
        # E progresses to I
        progress = Change(EXPOSED, INFECTIOUS)
        # I recovers to R
        recover = Change(INFECTIOUS, RECOVERED)
        move_s = Move(SUSCEPTIBLE)
        move_e = Move(EXPOSED)
        move_i = Move(INFECTIOUS)
        move_r = Move(RECOVERED)
        return [infect, progress, recover, move_s, move_e, move_i, move_r]

    def _build_network(self, params):
        raise NotImplementedError

    def _get_patch_seeding(self, params):
        seed = {}
        for n in self.network().nodes():
            seed[n] = {TypedMetapopulationNetwork.COMPARTMENTS: {SUSCEPTIBLE: params[SEIRDynamics.INIT_S],
                                                                 INFECTIOUS: params[SEIRDynamics.INIT_I]}}
        return seed

    def _get_edge_seeding(self, params):
        return {}