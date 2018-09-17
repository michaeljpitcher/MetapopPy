import unittest
from metapoppy import *

compartments = ['a','b','c']
patch_attributes = ['d','e','f']
edge_attributes = ['g','h','i']


class NAEvent1(Event):
    def __init__(self):
        Event.__init__(self)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, compartments[0])

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {compartments[1]: 1})


class NAEvent2(Event):
    def __init__(self):
        Event.__init__(self)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, compartments[1]) * 2

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {compartments[2]: 2})

class NADynamics(Dynamics):

    INITIAL_COMP_0 = 'initial_comp_0'
    INITIAL_COMP_1 = 'initial_comp_1'
    INITIAL_EDGE_0 = 'initial_edge_0'

    def __init__(self, g):
        Dynamics.__init__(self, g)

    def _create_events(self):
        events = [NAEvent1(), NAEvent2()]
        return events

    def _seed_network(self, params):
        value_comp_0 = params[NADynamics.INITIAL_COMP_0]
        value_comp_1 = params[NADynamics.INITIAL_COMP_1]
        for p in self._network.nodes():
            self._network.update_patch(p, {compartments[0]: value_comp_0, compartments[1]: value_comp_1})

    def _seed_events(self, event_parameters):
        event1 = next(e for e in self._events if isinstance(e, NAEvent1))
        event1.set_reaction_parameter(event_parameters[NAEvent1.__name__])
        event2 = next(e for e in self._events if isinstance(e, NAEvent2))
        event2.set_reaction_parameter(event_parameters[NAEvent2.__name__])

    def get_results(self):
        pass


class DynamicsTestCase(unittest.TestCase):

    def setUp(self):
        self.network = Network(compartments, patch_attributes, edge_attributes)
        self.network.add_nodes_from([1,2,3])
        self.network.add_edges_from([(1,2),(2,3)])

        self.dynamics = NADynamics(self.network)

    def test_setUp(self):
        params = {NAEvent1.__name__: 0.1, NAEvent2.__name__: 0.2, NADynamics.INITIAL_COMP_0: 3,
                  NADynamics.INITIAL_COMP_1: 5}
        self.dynamics.setUp(params)

        # Events set up
        rps = [0.1, 0.2]
        event1 = next(e for e in self.dynamics._events if isinstance(e, NAEvent1))
        self.assertEqual(event1._reaction_parameter, rps[0])
        event2 = next(e for e in self.dynamics._events if isinstance(e, NAEvent2))
        self.assertEqual(event2._reaction_parameter, rps[1])

        # Populations updated
        for a in [1, 2, 3]:
            self.assertEqual(self.dynamics._network.node[a][Network.COMPARTMENTS][compartments[0]],
                             params[NADynamics.INITIAL_COMP_0])
            self.assertEqual(self.dynamics._network.node[a][Network.COMPARTMENTS][compartments[1]],
                             params[NADynamics.INITIAL_COMP_1])
            self.assertEqual(self.dynamics._network.node[a][Network.COMPARTMENTS][compartments[2]], 0)

        # Population updates feeds through to event rates
        for col in range(0, 3):
            self.assertAlmostEqual(self.dynamics._rate_table[0][col], rps[0] * params[NADynamics.INITIAL_COMP_0])
            self.assertAlmostEqual(self.dynamics._rate_table[1][col], rps[1] * 2 * params[NADynamics.INITIAL_COMP_1])


if __name__ == '__main__':
    unittest.main()
