from epidemic import *
from networkx import *
from compartments import *
from ..events import *


class HerreraModel(Epidemic):
    """
    Herrera M, Bosch P, Najera M, Aguilera X. Modeling the spread of tuberculosis in semiclosed communities.
    Comput Math Methods Med 2013; 2013. DOI:10.1155/2013/648291.
    """

    INIT_S = 'initial_S'
    INIT_I = 'initial_II'

    def __init__(self):
        template = networkx.Graph()
        template.add_node(0)
        Epidemic.__init__(self, [SUSCEPTIBLE, EXPOSED, INFECTIOUS, NONINFECTIOUS, RECOVERED], template)

    def _create_events(self):
        # Birth
        birth = Create(SUSCEPTIBLE)
        # death S
        death_s = Death(SUSCEPTIBLE)
        # S to I_I
        s_to_ii = Infect(SUSCEPTIBLE, INFECTIOUS, INFECTIOUS)
        # S to I_N
        s_to_in = Infect(SUSCEPTIBLE, INFECTIOUS, NONINFECTIOUS)
        # S to E
        s_to_e = Infect(SUSCEPTIBLE, INFECTIOUS, EXPOSED)
        # death E
        death_e = Death(EXPOSED)
        # E to II delta
        e_to_ii_delta = Infect(EXPOSED, INFECTIOUS, INFECTIOUS)
        # E to II nu
        e_to_ii_nu = Change(EXPOSED, INFECTIOUS)
        # E to IN delta
        e_to_in_delta = Infect(EXPOSED, INFECTIOUS, NONINFECTIOUS)
        # E to IN nu
        e_to_in_nu = Change(EXPOSED, NONINFECTIOUS)
        # mu death II
        death_ii = Death(INFECTIOUS)
        # II to R
        ii_to_r = Change(INFECTIOUS, RECOVERED)
        # mu death IN
        death_in = Death(NONINFECTIOUS)
        # IN to R c
        in_to_r = Change(NONINFECTIOUS, RECOVERED)
        # death R
        death_r = Death(RECOVERED)
        # R to II
        r_to_ii = Change(RECOVERED, INFECTIOUS)
        # R to IN
        r_to_in = Change(RECOVERED, NONINFECTIOUS)
        # R to E
        r_to_e = Infect(RECOVERED, INFECTIOUS, EXPOSED)

        return [birth, death_s, s_to_ii, s_to_in, s_to_e,
                e_to_ii_delta, e_to_ii_nu, e_to_in_delta, e_to_in_nu, death_e,
                death_ii, ii_to_r,
                death_in, in_to_r,
                death_r, r_to_ii, r_to_in, r_to_e]

    def _build_network(self, params):
        raise NotImplementedError

    def _get_initial_patch_seeding(self, params):
        seed = {0: {TypedEnvironment.COMPARTMENTS: {SUSCEPTIBLE: params[HerreraModel.INIT_S],
                                                    INFECTIOUS: params[HerreraModel.INIT_I]}}}
        return seed

    def _seed_activated_patch(self, patch_id, params):
        seed = {}
        return seed

    def _get_initial_edge_seeding(self, params):
        return {}
