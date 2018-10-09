import unittest
from metapoppy import *

compartments = ['a','b','c']
patch_attributes = ['d','e','f']
edge_attributes = ['g','h','i']

RP_1_KEY = 'rp1'
RP_2_KEY = 'rp2'


class NAEvent1(Event):
    def __init__(self, rp_key):
        Event.__init__(self, rp_key)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, compartments[0])

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {compartments[1]: 1})


class NAEvent2(Event):
    def __init__(self, rp_key):
        Event.__init__(self, rp_key)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, compartments[1]) * 2

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {compartments[2]: 2})


class NADynamics(Dynamics):

    INITIAL_COMP_0 = 'initial_comp_0'
    INITIAL_COMP_1 = 'initial_comp_1'
    INITIAL_EDGE_0 = 'initial_edge_0'

    def __init__(self, network):
        Dynamics.__init__(self, network)

    def _create_events(self):
        events = [NAEvent1(RP_1_KEY), NAEvent2(RP_2_KEY)]
        return events

    def _seed_prototype_network(self, params):
        value_comp_0 = params[NADynamics.INITIAL_COMP_0]
        value_comp_1 = params[NADynamics.INITIAL_COMP_1]
        for p in self._network_prototype.nodes():
            self._network_prototype.update_patch(p, {compartments[0]: value_comp_0, compartments[1]: value_comp_1})

    def _seed_events(self, params):
        event1 = next(e for e in self._events if isinstance(e, NAEvent1))
        event1.set_parameters(params)
        event2 = next(e for e in self._events if isinstance(e, NAEvent2))
        event2.set_parameters(params)

    def _get_results(self):
        pass


class DynamicsTestCase(unittest.TestCase):

    def setUp(self):
        self.network = Network(compartments, patch_attributes, edge_attributes)
        self.network.add_nodes_from(['a','b','c'])
        self.network.add_edges_from([('a','b'),('b','c')])

        self.dynamics = NADynamics(self.network)

    def test_configure(self):
        params = {RP_1_KEY: 0.1, RP_2_KEY: 0.2, NADynamics.INITIAL_COMP_0: 3, NADynamics.INITIAL_COMP_1: 5}
        self.dynamics.configure(params)

        # Events set up
        event1 = next(e for e in self.dynamics._events if isinstance(e, NAEvent1))
        self.assertEqual(event1._reaction_parameter_key, RP_1_KEY)
        self.assertEqual(event1._reaction_parameter, params[RP_1_KEY])
        event2 = next(e for e in self.dynamics._events if isinstance(e, NAEvent2))
        self.assertEqual(event2._reaction_parameter_key, RP_2_KEY)
        self.assertEqual(event2._reaction_parameter, params[RP_2_KEY])

    def test_setUp(self):
        params = {RP_1_KEY: 0.1, RP_2_KEY: 0.2, NADynamics.INITIAL_COMP_0: 3,
                  NADynamics.INITIAL_COMP_1: 5}
        self.dynamics.configure(params)
        self.dynamics.setUp(params)

        # Populations updated
        for a in ['a','b','c']:
            self.assertEqual(self.dynamics._network.node[a][Network.COMPARTMENTS][compartments[0]],
                             params[NADynamics.INITIAL_COMP_0])
            self.assertEqual(self.dynamics._network.node[a][Network.COMPARTMENTS][compartments[1]],
                             params[NADynamics.INITIAL_COMP_1])
            self.assertEqual(self.dynamics._network.node[a][Network.COMPARTMENTS][compartments[2]], 0)

        # Population updates feeds through to event rates
        for col in range(0, 3):
            self.assertAlmostEqual(self.dynamics._rate_table[0][col], params[RP_1_KEY] *
                                   params[NADynamics.INITIAL_COMP_0])
            self.assertAlmostEqual(self.dynamics._rate_table[1][col], params[RP_2_KEY] * 2 *
                                   params[NADynamics.INITIAL_COMP_1])

    def test_run(self):
        params = {RP_1_KEY: 0.1, RP_2_KEY: 0.2, NADynamics.INITIAL_COMP_0: 3,
                  NADynamics.INITIAL_COMP_1: 5}
        self.dynamics.configure(params)
        self.dynamics.setUp(params)


if __name__ == '__main__':
    unittest.main()
