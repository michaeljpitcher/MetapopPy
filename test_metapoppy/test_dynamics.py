import unittest
from metapoppy import *

compartments = ['a','b','c']
patch_attributes = ['d','e','f']
edge_attributes = ['g','h','i']


class NAEvent1(Event):
    RP_1_KEY = 'rp1'
    def __init__(self):
        Event.__init__(self, [compartments[0]], [], [])

    def _define_parameter_keys(self):
        return NAEvent1.RP_1_KEY, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, compartments[0])

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {compartments[1]: 1})


class NAEvent2(Event):
    RP_2_KEY = 'rp2'
    def __init__(self):
        Event.__init__(self, [compartments[1]], [], [])

    def _define_parameter_keys(self):
        return NAEvent2.RP_2_KEY, []

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
        self.event1 = NAEvent1()
        self.event2 = NAEvent2()
        return [self.event1, self.event2]

    def _get_initial_patch_seeding(self, params):
        seeding = {}
        value_comp_0 = params[NADynamics.INITIAL_COMP_0]
        value_comp_1 = params[NADynamics.INITIAL_COMP_1]
        value_att_0 = params[NADynamics.INITIAL_ATT_0]
        value_att_1 = params[NADynamics.INITIAL_ATT_1]

        for p in self._network.nodes():
            seeding[p] = {Environment.COMPARTMENTS: {compartments[0]: value_comp_0, compartments[1]: value_comp_1},
                          Environment.ATTRIBUTES: {patch_attributes[0]: value_att_0, patch_attributes[1]: value_att_1}}

        return seeding

    def _seed_activated_patch(self, patch_id, params):
        return {}

    def _get_initial_edge_seeding(self, params):
        seeding = {}
        for u,v in self._network.edges():
            seeding[(u,v)] = {edge_attributes[0]: params[NADynamics.INITIAL_EDGE_0],
                              edge_attributes[1]: params[NADynamics.INITIAL_EDGE_1]}
        return seeding


class DynamicsTestCase(unittest.TestCase):

    def setUp(self):
        self.network = Environment(compartments, patch_attributes, edge_attributes)
        self.nodes = ['a1','b1','c1']
        self.network.add_nodes_from(self.nodes)
        self.network.add_edges_from([('a1','b1'),('b1','c1')])

        self.dynamics = NADynamics(self.network)

    def test_required_event_parameters(self):
        self.assertItemsEqual(self.dynamics.required_event_parameters(), [NAEvent1.RP_1_KEY, NAEvent2.RP_2_KEY])

    def test_set_start_time_maximum_time_record_interval(self):

        params = {NAEvent1.RP_1_KEY: 0.1, NAEvent2.RP_2_KEY: 0.2, NADynamics.INITIAL_COMP_0: 3,
                  NADynamics.INITIAL_COMP_1: 5,
                  NADynamics.INITIAL_ATT_0: 7, NADynamics.INITIAL_ATT_1: 11, NADynamics.INITIAL_EDGE_0: 13,
                  NADynamics.INITIAL_EDGE_1: 17}
        self.dynamics.set_start_time(101.0)
        self.dynamics.set_maximum_time(102.0)
        self.dynamics.set_record_interval(0.25)

        self.dynamics.configure(params)
        self.dynamics.setUp(params)
        res = self.dynamics.do(params)

        self.assertItemsEqual(res.keys(), [101.0, 101.25, 101.5, 101.75, 102.0])

    def test_configure(self):
        params = {NAEvent1.RP_1_KEY: 0.1, NAEvent2.RP_2_KEY: 0.2, NADynamics.INITIAL_COMP_0: 3,
                  NADynamics.INITIAL_COMP_1: 5,
                  NADynamics.INITIAL_ATT_0: 7, NADynamics.INITIAL_ATT_1: 11, NADynamics.INITIAL_EDGE_0: 13,
                  NADynamics.INITIAL_EDGE_1: 17}
        self.dynamics.configure(params)

        # Events set up
        event1 = next(e for e in self.dynamics._events if isinstance(e, NAEvent1))
        self.assertEqual(event1._reaction_parameter_key, NAEvent1.RP_1_KEY)
        self.assertEqual(event1._reaction_parameter, params[NAEvent1.RP_1_KEY])
        event2 = next(e for e in self.dynamics._events if isinstance(e, NAEvent2))
        self.assertEqual(event2._reaction_parameter_key, NAEvent2.RP_2_KEY)
        self.assertEqual(event2._reaction_parameter, params[NAEvent2.RP_2_KEY])

    def test_setUp(self):
        params = {NAEvent1.RP_1_KEY: 0.1, NAEvent2.RP_2_KEY: 0.2, NADynamics.INITIAL_COMP_0: 3, NADynamics.INITIAL_COMP_1: 5,
                  NADynamics.INITIAL_ATT_0: 7, NADynamics.INITIAL_ATT_1: 11, NADynamics.INITIAL_EDGE_0: 13,
                  NADynamics.INITIAL_EDGE_1: 17}
        self.dynamics.configure(params)
        self.dynamics.setUp(params)

        # Populations updated - ensure seeding from configure used in setup
        for a in ['a1','b1','c1']:
            self.assertEqual(self.dynamics._network.node[a][Environment.COMPARTMENTS][compartments[0]],
                             params[NADynamics.INITIAL_COMP_0])
            self.assertEqual(self.dynamics._network.node[a][Environment.COMPARTMENTS][compartments[1]],
                             params[NADynamics.INITIAL_COMP_1])
            self.assertEqual(self.dynamics._network.node[a][Environment.COMPARTMENTS][compartments[2]], 0)

        # Population updates feeds through to event rates
        # for col in range(0, 3):
        #     self.assertAlmostEqual(self.dynamics._rate_table[col][0], params[NAEvent1.RP_1_KEY] *
        #                            params[NADynamics.INITIAL_COMP_0])
        #     self.assertAlmostEqual(self.dynamics._rate_table[col][1], params[NAEvent2.RP_2_KEY] * 2 *
        #                            params[NADynamics.INITIAL_COMP_1])

        # Attribute updates
        for a in ['a1','b1','c1']:
            self.assertEqual(self.dynamics._network.node[a][Environment.ATTRIBUTES][patch_attributes[0]],
                             params[NADynamics.INITIAL_ATT_0])
            self.assertEqual(self.dynamics._network.node[a][Environment.ATTRIBUTES][patch_attributes[1]],
                             params[NADynamics.INITIAL_ATT_1])
            self.assertEqual(self.dynamics._network.node[a][Environment.ATTRIBUTES][patch_attributes[2]], 0)

        # Edge updates
        for _, _, d in self.dynamics._network.edges(data=True):
            self.assertEqual(d[edge_attributes[0]], params[NADynamics.INITIAL_EDGE_0])
            self.assertEqual(d[edge_attributes[1]], params[NADynamics.INITIAL_EDGE_1])
            self.assertFalse(d[edge_attributes[2]])

    def test_run(self):
        params = {NAEvent1.RP_1_KEY: 0.1, NAEvent2.RP_2_KEY: 0.2, NADynamics.INITIAL_COMP_0: 3, NADynamics.INITIAL_COMP_1: 5,
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
                self.assertItemsEqual(q.keys(), [Environment.ATTRIBUTES, Environment.COMPARTMENTS])
                self.assertItemsEqual(q[Environment.COMPARTMENTS], compartments)
                self.assertItemsEqual(q[Environment.ATTRIBUTES], patch_attributes)

    def test_do(self):
        params = {NAEvent1.RP_1_KEY: 0.1, NAEvent2.RP_2_KEY: 0.2,
                  NADynamics.INITIAL_COMP_0: 3, NADynamics.INITIAL_COMP_1: 5,
                  NADynamics.INITIAL_ATT_0: 7, NADynamics.INITIAL_ATT_1: 11,
                  NADynamics.INITIAL_EDGE_0: 13, NADynamics.INITIAL_EDGE_1: 17}
        self.dynamics.configure(params)
        self.dynamics.setUp(params)
        self.dynamics.do(params)

    def test_post_event(self):
        letters = ['a','b','c','d','e','f','g','h','i']
        strq = ''.join(numpy.random.choice(letters, 10, True))

        def func(atts):
            return strq

        self.dynamics.post_event(14.5, lambda a: func(a), [])
        self.assertEqual(len(self.dynamics._posted_events), 1)
        self.assertEqual(self.dynamics._posted_events[0][0], 14.5)

        # Check we get the random string back
        self.assertEqual(self.dynamics._posted_events[0][1](None), strq)


class EventNoDep(Event):
    def __init__(self):
        Event.__init__(self, [], [], [])

    def _define_parameter_keys(self):
        return self.__class__.__name__, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return 1

    def perform(self, network, patch_id):
        pass


class EventPatchCompDep(Event):
    def __init__(self, comp):
        self.comp = comp
        Event.__init__(self, [comp], [], [])

    def _define_parameter_keys(self):
        return self.__class__.__name__, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self.comp)

    def perform(self, network, patch_id):
        pass


class EventPatchAttDep(Event):
    def __init__(self, att):
        self.att = att
        Event.__init__(self, [], [att], [])

    def _define_parameter_keys(self):
        return self.__class__.__name__, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_attribute_value(patch_id, self.att)

    def perform(self, network, patch_id):
        pass


class EventEdgeAttDep(Event):
    def __init__(self, att):
        self.att = att
        Event.__init__(self, [], [], [att])

    def _define_parameter_keys(self):
        return self.__class__.__name__, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return sum(v[self.att] for v in network[patch_id].values())

    def perform(self, network, patch_id):
        pass


class PropDynamics(Dynamics):
    def __init__(self, network):
        Dynamics.__init__(self, network)

    def _create_events(self):
        events = [EventNoDep(), EventPatchCompDep(compartments[0]), EventPatchAttDep(patch_attributes[0]),
                  EventEdgeAttDep(edge_attributes[0])]
        return events

    def _get_initial_patch_seeding(self, params):
        return {}

    def _get_initial_edge_seeding(self, params):
        return {}

    def _seed_activated_patch(self, patch_id, params):
        return {}


class PropagationTestCase(unittest.TestCase):

    def setUp(self):
        self.network = Environment(compartments, patch_attributes, edge_attributes)
        self.nodes = ['a1', 'b1', 'c1']
        self.network.add_nodes_from(self.nodes)
        self.network.add_edges_from([('a1', 'b1'), ('b1', 'c1')])
        self.dynamics = PropDynamics(self.network)
#
    def test_propagation(self):
        params = {EventNoDep.__name__: 0.1, EventPatchCompDep.__name__: 0.2,
                  EventPatchAttDep.__name__: 0.3, EventEdgeAttDep.__name__: 0.4}
        self.dynamics.configure(params)
        self.dynamics.setUp(params)

        # No values so check rates - NoDep should have value, rest should be 0
        nodep_rates = self.dynamics._rate_table[:,0]
        patchcompdep_rates = self.dynamics._rate_table[:,1]
        patchattdep_rates = self.dynamics._rate_table[:,2]
        edgeattdep_rates = self.dynamics._rate_table[:,3]
        for r in nodep_rates:
            self.assertEqual(r, params[EventNoDep.__name__] * 1)
        for r in patchcompdep_rates:
            self.assertFalse(r)
        for r in patchattdep_rates:
            self.assertFalse(r)
        for r in edgeattdep_rates:
            self.assertFalse(r)

        # Update compartment
        self.network.update_patch(self.nodes[0], {compartments[0]: 1})
        self.network.update_patch(self.nodes[1], {compartments[0]: 2})
        self.network.update_patch(self.nodes[2], {compartments[0]: 3})

        nodep_rates = self.dynamics._rate_table[:,0]
        patchcompdep_rates = self.dynamics._rate_table[:,1]
        patchattdep_rates = self.dynamics._rate_table[:,2]
        edgeattdep_rates = self.dynamics._rate_table[:,3]
        for r in nodep_rates:
            self.assertEqual(r, params[EventNoDep.__name__] * 1)
        for i in range(len(patchcompdep_rates)):
            r = patchcompdep_rates[i]
            p = self.dynamics._active_patches[i]
            self.assertEqual(r, params[EventPatchCompDep.__name__] * self.network.node[p][Environment.COMPARTMENTS][compartments[0]])
        for r in patchattdep_rates:
            self.assertFalse(r)
        for r in edgeattdep_rates:
            self.assertFalse(r)

        # Update attribute
        self.network.update_patch(self.nodes[0], attribute_changes={patch_attributes[0]: 4})
        self.network.update_patch(self.nodes[1], attribute_changes={patch_attributes[0]: 5})
        self.network.update_patch(self.nodes[2], attribute_changes={patch_attributes[0]: 6})

        nodep_rates = self.dynamics._rate_table[:,0]
        patchcompdep_rates = self.dynamics._rate_table[:,1]
        patchattdep_rates = self.dynamics._rate_table[:,2]
        edgeattdep_rates = self.dynamics._rate_table[:,3]
        for r in nodep_rates:
            self.assertEqual(r, params[EventNoDep.__name__] * 1)
        for i in range(len(patchcompdep_rates)):
            r = patchcompdep_rates[i]
            p = self.dynamics._active_patches[i]
            self.assertEqual(r, params[EventPatchCompDep.__name__] * self.network.node[p][Environment.COMPARTMENTS][compartments[0]])
        for i in range(len(patchattdep_rates)):
            r = patchattdep_rates[i]
            p = self.dynamics._active_patches[i]
            self.assertEqual(r, params[EventPatchAttDep.__name__] * self.network.node[p][Environment.ATTRIBUTES][patch_attributes[0]])
        for r in edgeattdep_rates:
            self.assertFalse(r)

        self.network.update_edge(self.nodes[0], self.nodes[1], {edge_attributes[0]: 11})
        self.network.update_edge(self.nodes[1], self.nodes[2], {edge_attributes[0]: 13})

        for r in nodep_rates:
            self.assertEqual(r, params[EventNoDep.__name__] * 1)
        for i in range(len(patchcompdep_rates)):
            r = patchcompdep_rates[i]
            p = self.dynamics._active_patches[i]
            self.assertEqual(r, params[EventPatchCompDep.__name__] * self.network.node[p][Environment.COMPARTMENTS][compartments[0]])
        for i in range(len(patchattdep_rates)):
            r = patchattdep_rates[i]
            p = self.dynamics._active_patches[i]
            self.assertEqual(r, params[EventPatchAttDep.__name__] * self.network.node[p][Environment.ATTRIBUTES][patch_attributes[0]])
        for i in range(len(edgeattdep_rates)):
            r = edgeattdep_rates[i]
            p = self.dynamics._active_patches[i]
            v = sum([a[edge_attributes[0]] for a in self.network[p].values()])
            self.assertEqual(r, params[EventEdgeAttDep.__name__] * v)


if __name__ == '__main__':
    unittest.main()
