from metapoppy import *
from epidemic import *
from ..events import *
from compartments import *


class SISDynamics(Epidemic):

    RP_INFECT = 'infection_rate'
    RP_RECOVER = 'recover_rate'
    RP_MOVE_S = 's_move_rate'
    RP_MOVE_I = 'i_move_rate'

    INIT_S = 'initial_population_susceptible'
    INIT_I = 'initial_population_infected'

    def __init__(self, template_network):
        Epidemic.__init__(self, [SUSCEPTIBLE, INFECTIOUS], template_network)

    def _create_events(self):
        infect = Infect(SUSCEPTIBLE, INFECTIOUS, INFECTIOUS)
        recover = Change(INFECTIOUS, SUSCEPTIBLE)
        move_s = Move(SUSCEPTIBLE)
        move_i = Move(INFECTIOUS)
        return [infect, recover, move_s, move_i]

    def _build_network(self, params):
        raise NotImplementedError

    def _get_patch_seeding(self, params):
        seed = {}
        for n in self.network().nodes():
            seed[n] = {TypedMetapopulationNetwork.COMPARTMENTS: {SUSCEPTIBLE: params[SISDynamics.INIT_S],
                                                                 INFECTIOUS: params[SISDynamics.INIT_I]}}
        return seed

    def _get_edge_seeding(self, params):
        return {}