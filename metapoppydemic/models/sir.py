from metapoppy import *
from epidemic import *
from ..events import *
from compartments import *


class SIRDynamics(Epidemic):

    INIT_S = 'initial_population_susceptible'
    INIT_I = 'initial_population_infected'

    def __init__(self, template_network):
        Epidemic.__init__(self, [SUSCEPTIBLE, INFECTIOUS, RECOVERED], template_network)

    def _create_events(self):
        infect = Infect(SUSCEPTIBLE, INFECTIOUS, INFECTIOUS)
        recover = Change(INFECTIOUS, RECOVERED)
        move_s = Move(SUSCEPTIBLE)
        move_i = Move(INFECTIOUS)
        move_r = Move(RECOVERED)
        return [infect, recover, move_s, move_i, move_r]

    def _build_network(self, params):
        raise NotImplementedError

    def _get_patch_seeding(self, params):
        seed = {}
        for n in self.network().nodes():
            seed[n] = {TypedMetapopulationNetwork.COMPARTMENTS: {SUSCEPTIBLE: params[SIRDynamics.INIT_S],
                                                                 INFECTIOUS: params[SIRDynamics.INIT_I]}}
        return seed

    def _get_edge_seeding(self, params):
        return {}