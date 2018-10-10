from metapoppy import *
from epidemic import *
from ..events import *
from compartments import *


class SEIRDynamics(Epidemic):

    RP_INFECT = 'infection_rate'
    RP_PROGRESS = 'progress_rate'
    RP_RECOVER = 'recover_rate'
    RP_MOVE_S = 's_move_rate'
    RP_MOVE_E = 'e_move_rate'
    RP_MOVE_I = 'i_move_rate'
    RP_MOVE_R = 'r_move_rate'

    INIT_S = 'initial_population_susceptible'
    INIT_I = 'initial_population_infected'

    def __init__(self, template_network):
        Epidemic.__init__(self, [SUSCEPTIBLE, EXPOSED, INFECTIOUS, RECOVERED], template_network)

    def _create_events(self):
        # Contact with I moves S to E
        infect = Infect(SEIRDynamics.RP_INFECT, SUSCEPTIBLE, INFECTIOUS, EXPOSED)
        # E progresses to I
        progress = Change(SEIRDynamics.RP_PROGRESS, EXPOSED, INFECTIOUS)
        # I recovers to R
        recover = Change(SEIRDynamics.RP_RECOVER, INFECTIOUS, RECOVERED)
        move_s = Move(SEIRDynamics.RP_MOVE_S, SUSCEPTIBLE)
        move_e = Move(SEIRDynamics.RP_MOVE_E, EXPOSED)
        move_i = Move(SEIRDynamics.RP_MOVE_I, INFECTIOUS)
        move_r = Move(SEIRDynamics.RP_MOVE_R, RECOVERED)
        return [infect, progress, recover, move_s, move_e, move_i, move_r]

    def _seed_prototype_network(self, params):
        for n in self.network().nodes():
            self._network.update_patch(n, {SUSCEPTIBLE: params[SEIRDynamics.INIT_S]})
            self._network.update_patch(n, {INFECTIOUS: params[SEIRDynamics.INIT_I]})
