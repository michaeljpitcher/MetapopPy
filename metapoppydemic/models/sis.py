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
        infect = Infect(SISDynamics.RP_INFECT, SUSCEPTIBLE, INFECTIOUS, INFECTIOUS)
        recover = Change(SISDynamics.RP_RECOVER, INFECTIOUS, SUSCEPTIBLE)
        move_s = Move(SISDynamics.RP_MOVE_S, SUSCEPTIBLE)
        move_i = Move(SISDynamics.RP_MOVE_I, INFECTIOUS)
        return [infect, recover, move_s, move_i]

    def _seed_prototype_network(self, params):
        for n in self.network().nodes():
            self._network.update_patch(n, {SUSCEPTIBLE: params[SISDynamics.INIT_S]})
            self._network.update_patch(n, {INFECTIOUS: params[SISDynamics.INIT_I]})
