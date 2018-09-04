import unittest
from metapoppy_group_by_patch import *


class NAEvent(Event):
    REACTION_PARAMETER = 0.5

    def __init__(self, sv_comp, perf_comp):
        self.sv_comp = sv_comp
        self.perf_comp = perf_comp
        Event.__init__(self)

    def _calculate_state_variable(self, patch, edges):
        return patch[Network.COMPARTMENTS][self.sv_comp]

    def perform(self, network):
        network.update_patch(self._patch_id, {self.perf_comp: 2})


class NADynamics(Dynamics):
    def __init__(self, g):
        Dynamics.__init__(self, g)

    def _create_events(self, network):
        event1 = NAEvent(self._compartments[0], self._compartments[1])
        network.attach_event(event1, 1)
        event2 = NAEvent(self._compartments[1], self._compartments[2])
        network.attach_event(event2, 2)

    def _seed_events(self, params):
        NAEvent.REACTION_PARAMETER = params[NAEvent.__name__][Event.REACTION_PARAMETER]


class DynamicsTestCase(unittest.TestCase):

    def setUp(self):
        self.comps = ['a','b','c']
        self.p_atts = ['d','e','f']
        self.e_atts = ['g','h','i']

        self.network = Network(self.comps, self.p_atts, self.e_atts)
        self.network.add_nodes_from([1,2,3])
        self.network.add_edges_from([(1,2),(2,3)])

        self.dynamics = NADynamics(self.network)

    def test_setUp(self):
        params = {Dynamics.INITIAL_PATCHES: {1: {Network.COMPARTMENTS:{self.comps[0]: 1},
                                                 Network.ATTRIBUTES:{self.p_atts[0]: 0.5, self.p_atts[1]: 1.5,
                                                                     self.p_atts[2]: 3.5}},
                                             2: {Network.COMPARTMENTS:{self.comps[1]: 6},
                                                 Network.ATTRIBUTES:{self.p_atts[0]: 0.5, self.p_atts[1]: 1.5,
                                                                     self.p_atts[2]: 3.5}}},
                  Dynamics.INITIAL_EDGES: {(1, 2): {self.e_atts[0]: 0.1, self.e_atts[1]: 0.2, self.e_atts[2]: 0.3},
                                           (2, 3): {self.e_atts[0]: 0.4, self.e_atts[1]: 0.5, self.e_atts[2]: 0.6}},
                  Network.EVENTS: {NAEvent.__name__: {Event.REACTION_PARAMETER: 0.11}}
                  }
        self.dynamics.setUp(params)

        # Graph
        self.assertTrue(self.dynamics.network())

        # Events
        self.assertEqual(len(self.dynamics.network().node[1][Network.EVENTS]), 1)
        event_at_1 = self.dynamics.network().node[1][Network.EVENTS][0]
        self.assertTrue(isinstance(event_at_1, NAEvent))
        self.assertEqual(event_at_1.sv_comp, self.comps[0])
        self.assertEqual(event_at_1.perf_comp, self.comps[1])
        self.assertEqual(event_at_1.REACTION_PARAMETER, 0.11)

        self.assertEqual(len(self.dynamics.network().node[2][Network.EVENTS]), 1)
        event_at_2 = self.dynamics.network().node[2][Network.EVENTS][0]
        self.assertTrue(isinstance(event_at_2, NAEvent))
        self.assertEqual(event_at_2.sv_comp, self.comps[1])
        self.assertEqual(event_at_2.perf_comp, self.comps[2])
        self.assertEqual(event_at_2.REACTION_PARAMETER, 0.11)

        self.assertEqual(len(self.dynamics.network().node[3][Network.EVENTS]), 0)

        # Comps & atts
        for p, data in params[Dynamics.INITIAL_PATCHES].iteritems():
            for c in self.comps:
                if c in data[Network.COMPARTMENTS]:
                    self.assertEqual(self.dynamics.network().get_compartment_value(p, c), data[Network.COMPARTMENTS][c])
                else:
                    self.assertEqual(self.dynamics.network().get_compartment_value(p, c), 0)
            for a in self.p_atts:
                if a in data[Network.ATTRIBUTES]:
                    self.assertEqual(self.dynamics.network().get_attribute_value(p, a), data[Network.ATTRIBUTES][a])
                else:
                    self.assertEqual(self.dynamics.network().get_attribute_value(p, a), 0)

        # Edges
        for (u, v), data in params[Dynamics.INITIAL_EDGES].iteritems():
            edge = self.dynamics.network().edge[u][v]
            for a in self.e_atts:
                self.assertEqual(edge[a], data[a])


if __name__ == '__main__':
    unittest.main()
