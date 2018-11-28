import unittest
from metapoppy import *

compartments = ['a','b','c']
patch_attributes = ['d','e','f']
edge_attributes = ['g','h','i']

RP_1_KEY = 'rp1'
RP_2_KEY = 'rp2'


class NAEvent1(Event):
    def __init__(self, rp_key):
        Event.__init__(self, [compartments[0]], [], rp_key)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, compartments[0])

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {compartments[1]: 1})


class NAEvent2(Event):
    def __init__(self, rp_key):
        Event.__init__(self, [compartments[1]], [], rp_key)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, compartments[1]) * 2

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {compartments[2]: 2})


class NADynamics(Dynamics):

    INITIAL_COMP_0 = 'initial_comp_0'
    INITIAL_COMP_1 = 'initial_comp_1'
    INITIAL_ATT_0 = 'initial_att_0'
    INITIAL_ATT_1 = 'initial_att_1'
    INITIAL_EDGE_0 = 'initial_edge_0'
    INITIAL_EDGE_1 = 'initial_edge_1'

    def __init__(self, network):
        Dynamics.__init__(self, network)

    def _create_events(self):
        self.event1 = NAEvent1(RP_1_KEY)
        self.event2 = NAEvent2(RP_2_KEY)
        return [self.event1, self.event2]

    def _get_patch_seeding(self, params):
        seeding = {}
        value_comp_0 = params[NADynamics.INITIAL_COMP_0]
        value_comp_1 = params[NADynamics.INITIAL_COMP_1]
        value_att_0 = params[NADynamics.INITIAL_ATT_0]
        value_att_1 = params[NADynamics.INITIAL_ATT_1]

        for p in self._network.nodes():
            seeding[p] = {Network.COMPARTMENTS: {compartments[0]: value_comp_0, compartments[1]: value_comp_1},
                          Network.ATTRIBUTES: {patch_attributes[0]: value_att_0, patch_attributes[1]: value_att_1}}

        return seeding

    def _get_edge_seeding(self, params):
        seeding = {}
        for u,v in self._network.edges():
            seeding[(u,v)] = {edge_attributes[0]: params[NADynamics.INITIAL_EDGE_0],
                              edge_attributes[1]: params[NADynamics.INITIAL_EDGE_1]}
        return seeding


class DynamicsTestCase(unittest.TestCase):

    def setUp(self):
        self.network = Network(compartments, patch_attributes, edge_attributes)
        self.nodes = ['a1','b1','c1']
        self.network.add_nodes_from(self.nodes)
        self.network.add_edges_from([('a1','b1'),('b1','c1')])

        self.dynamics = NADynamics(self.network)

    def test_configure(self):
        params = {RP_1_KEY: 0.1, RP_2_KEY: 0.2, NADynamics.INITIAL_COMP_0: 3, NADynamics.INITIAL_COMP_1: 5,
                  NADynamics.INITIAL_ATT_0: 7, NADynamics.INITIAL_ATT_1: 11, NADynamics.INITIAL_EDGE_0: 13,
                  NADynamics.INITIAL_EDGE_1: 17}
        self.dynamics.configure(params)

        # Events set up
        event1 = next(e for e in self.dynamics._events if isinstance(e, NAEvent1))
        self.assertEqual(event1._reaction_parameter_key, RP_1_KEY)
        self.assertEqual(event1._reaction_parameter, params[RP_1_KEY])
        event2 = next(e for e in self.dynamics._events if isinstance(e, NAEvent2))
        self.assertEqual(event2._reaction_parameter_key, RP_2_KEY)
        self.assertEqual(event2._reaction_parameter, params[RP_2_KEY])

    def test_setUp(self):
        params = {RP_1_KEY: 0.1, RP_2_KEY: 0.2, NADynamics.INITIAL_COMP_0: 3, NADynamics.INITIAL_COMP_1: 5,
                  NADynamics.INITIAL_ATT_0: 7, NADynamics.INITIAL_ATT_1: 11, NADynamics.INITIAL_EDGE_0: 13,
                  NADynamics.INITIAL_EDGE_1: 17}
        self.dynamics.configure(params)
        self.dynamics.setUp(params)

        # Populations updated - ensure seeding from configure used in setup
        for a in ['a1','b1','c1']:
            self.assertEqual(self.dynamics._network.node[a][Network.COMPARTMENTS][compartments[0]],
                             params[NADynamics.INITIAL_COMP_0])
            self.assertEqual(self.dynamics._network.node[a][Network.COMPARTMENTS][compartments[1]],
                             params[NADynamics.INITIAL_COMP_1])
            self.assertEqual(self.dynamics._network.node[a][Network.COMPARTMENTS][compartments[2]], 0)

        # Population updates feeds through to event rates
        for col in range(0, 3):
            self.assertAlmostEqual(self.dynamics._rate_table[col][0], params[RP_1_KEY] *
                                   params[NADynamics.INITIAL_COMP_0])
            self.assertAlmostEqual(self.dynamics._rate_table[col][1], params[RP_2_KEY] * 2 *
                                   params[NADynamics.INITIAL_COMP_1])

        # Attribute updates
        for a in ['a1','b1','c1']:
            self.assertEqual(self.dynamics._network.node[a][Network.ATTRIBUTES][patch_attributes[0]],
                             params[NADynamics.INITIAL_ATT_0])
            self.assertEqual(self.dynamics._network.node[a][Network.ATTRIBUTES][patch_attributes[1]],
                             params[NADynamics.INITIAL_ATT_1])
            self.assertEqual(self.dynamics._network.node[a][Network.ATTRIBUTES][patch_attributes[2]], 0)

        # Edge updates
        for _, _, d in self.dynamics._network.edges(data=True):
            self.assertEqual(d[edge_attributes[0]], params[NADynamics.INITIAL_EDGE_0])
            self.assertEqual(d[edge_attributes[1]], params[NADynamics.INITIAL_EDGE_1])
            self.assertFalse(d[edge_attributes[2]])

    def test_run(self):
        params = {RP_1_KEY: 0.1, RP_2_KEY: 0.2, NADynamics.INITIAL_COMP_0: 3, NADynamics.INITIAL_COMP_1: 5,
                  NADynamics.INITIAL_ATT_0: 7, NADynamics.INITIAL_ATT_1: 11, NADynamics.INITIAL_EDGE_0: 13,
                  NADynamics.INITIAL_EDGE_1: 17}
        self.dynamics.set(params)

        max_time = 69
        interval = 0.1
        self.dynamics.set_maximum_time(max_time)
        self.dynamics.set_record_interval(interval)

        r = self.dynamics.run()
        self.assertItemsEqual(r.keys(), ['results', 'metadata', 'parameters'])
        self.assertItemsEqual(r['parameters'].keys(), params.keys())
        for k,v in r['parameters'].iteritems():
            self.assertEqual(v, params[k])
        self.assertEqual(len(r['results']), (max_time/interval)+1)
        for t,v in r['results'].iteritems():
            self.assertTrue(isinstance(t, float))
            self.assertItemsEqual(v.keys(), self.nodes)
            for n,q in v.iteritems():
                self.assertItemsEqual(q.keys(), [Network.ATTRIBUTES, Network.COMPARTMENTS])
                self.assertItemsEqual(q[Network.COMPARTMENTS], compartments)
                self.assertItemsEqual(q[Network.ATTRIBUTES], patch_attributes)

    def test_do(self):
        params = {RP_1_KEY: 0.1, RP_2_KEY: 0.2, NADynamics.INITIAL_COMP_0: 3, NADynamics.INITIAL_COMP_1: 5,
                  NADynamics.INITIAL_ATT_0: 7, NADynamics.INITIAL_ATT_1: 11, NADynamics.INITIAL_EDGE_0: 13,
                  NADynamics.INITIAL_EDGE_1: 17}
        self.dynamics.configure(params)
        self.dynamics.setUp(params)
        self.dynamics.do(params)


if __name__ == '__main__':
    unittest.main()
