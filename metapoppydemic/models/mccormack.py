from epidemic import *
from networkx import *
from compartments import *
from ..events import *
from ..environment import *


class McCormackModel(Dynamics):
    """
    McCormack RK, Allen LJS. Multi-patch deterministic and stochastic models for wildlife diseases.
    J Biol Dyn 2007; 1: 63-85.
    """

    INIT_S = 'initial_S'
    INIT_I = 'initial_I'
    INIT_R = 'initial_R'

    def __init__(self, b, d, lambdaa, k, alpha, gamma):
        self._b = b
        self._d = d
        self._lambdaa = lambdaa
        self._k = k
        self._alpha = alpha
        self._gamma = gamma
        env = McCormackEnvironment([SUSCEPTIBLE, INFECTIOUS, RECOVERED])
        Dynamics.__init__(self, env)

    def _create_events(self):
        all_comps = [SUSCEPTIBLE, INFECTIOUS, RECOVERED]

        # Birth
        birth = McCormackBirth(SUSCEPTIBLE, all_comps)
        # Death
        death_S = McCormackDeath(SUSCEPTIBLE, all_comps)
        death_I = McCormackDeath(INFECTIOUS, all_comps)
        death_R = McCormackDeath(RECOVERED, all_comps)
        # Infection
        infection = McCormackInfection(SUSCEPTIBLE, INFECTIOUS, all_comps)
        # Move
        move_s = McCormackMove(SUSCEPTIBLE)
        move_i = McCormackMove(INFECTIOUS)
        move_r = McCormackMove(RECOVERED)
        # Infection based death
        death_I_inf = McCormackDeathInfection(INFECTIOUS)
        # Recovery
        recover = McCormackRecover(INFECTIOUS, RECOVERED)

        return [birth, death_S, death_I, death_R, infection, move_s, move_i, move_r, death_I_inf, recover]

    def _build_network(self, params):
        raise NotImplementedError

    def _get_initial_patch_seeding(self, params):
        seed = {}
        for n in range(1,4):
            seed[n] = {TypedEnvironment.COMPARTMENTS: {SUSCEPTIBLE: params[McCormackModel.INIT_S],
                                                    INFECTIOUS: params[McCormackModel.INIT_I],
                                                    RECOVERED: params[McCormackModel.INIT_R]},
                       TypedEnvironment.ATTRIBUTES: {McCormackEnvironment.BIRTH_RATE: self._b[n],
                                                    McCormackEnvironment.BASE_DEATH_RATE:self._d[n][0],
                                                    McCormackEnvironment.POPULATION_DEATH_RATE:self._d[n][1],
                                                    McCormackEnvironment.INFECTION_LAMBDA:self._lambdaa[n],
                                                    McCormackEnvironment.CARRYING_CAPACITY:self._k[n],
                                                    McCormackEnvironment.INFECTION_DEATH_RATE:self._alpha[n],
                                                    McCormackEnvironment.RECOVERY_RATE:self._gamma[n]}}
        return seed

    def _seed_activated_patch(self, patch_id, params):
        return {}

    def _get_initial_edge_seeding(self, params):
        return {}