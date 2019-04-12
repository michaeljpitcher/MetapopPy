from metapoppy import *
from epidemic import *
from ..events import *
from compartments import *


class SISDynamics(Epidemic):

    INIT_S = 'initial_population_susceptible'
    INIT_I = 'initial_population_infected'

    def __init__(self, template_network):
        self.rp_infect_key = self.rp_recover_key = self.rp_move_s_key = self.rp_move_i_key = None
        Epidemic.__init__(self, [SUSCEPTIBLE, INFECTIOUS], template_network)

    def _create_events(self):
        infect = Infect(SUSCEPTIBLE, INFECTIOUS, INFECTIOUS)
        self.rp_infect_key = infect.reaction_parameter()
        recover = Change(INFECTIOUS, SUSCEPTIBLE)
        self.rp_recover_key = recover.reaction_parameter()
        move_s = Move(SUSCEPTIBLE)
        self.rp_move_s_key = move_s.reaction_parameter()
        move_i = Move(INFECTIOUS)
        self.rp_move_i_key = move_i.reaction_parameter()
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