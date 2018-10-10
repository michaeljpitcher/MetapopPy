from metapoppy import *
from epidemic import *
from ..events import *
from compartments import *


class SIRDynamics(Epidemic):

    RP_INFECT = 'infection_rate'
    RP_RECOVER = 'recover_rate'
    RP_MOVE_S = 's_move_rate'
    RP_MOVE_I = 'i_move_rate'
    RP_MOVE_R = 'r_move_rate'

    INIT_S = 'initial_population_susceptible'
    INIT_I = 'initial_population_infected'

    def __init__(self, template_network):
        Epidemic.__init__(self, [SUSCEPTIBLE, INFECTIOUS, RECOVERED], template_network)

    def _create_events(self):
        infect = Infect(SIRDynamics.RP_INFECT, SUSCEPTIBLE, INFECTIOUS, INFECTIOUS)
        recover = Change(SIRDynamics.RP_RECOVER, INFECTIOUS, RECOVERED)
        move_s = Move(SIRDynamics.RP_MOVE_S, SUSCEPTIBLE)
        move_i = Move(SIRDynamics.RP_MOVE_I, INFECTIOUS)
        move_r = Move(SIRDynamics.RP_MOVE_R, RECOVERED)
        return [infect, recover, move_s, move_i, move_r]

    def _seed_prototype_network(self, params):
        for n in self.network().nodes():
            self._network.update_patch(n, {SUSCEPTIBLE: params[SIRDynamics.INIT_S]})
            self._network.update_patch(n, {INFECTIOUS: params[SIRDynamics.INIT_I]})
