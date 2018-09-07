import unittest
from metapoppy_group_by_patch import *
import numpy


class NAEvent1(Event):
    def __init__(self, comp1, comp2):
        self.sv_comp = comp1
        self.perf_comp = comp2
        Event.__init__(self)

    def _calculate_state_variable(self, patch, edges):
        return patch[Network.COMPARTMENTS][self.sv_comp]

    def perform(self, network):
        network.update_patch(self._patch_id, {self.perf_comp: 1})


class NAEvent2(Event):
    def __init__(self, comp1, comp2):
        self.sv_comp = comp1
        self.perf_comp = comp2
        Event.__init__(self)

    def _calculate_state_variable(self, patch, edges):
        return 1.5 * patch[Network.COMPARTMENTS][self.sv_comp]

    def perform(self, network):
        network.update_patch(self._patch_id, {self.perf_comp: 1})


class NAEvent3(Event):
    def __init__(self, edge_att):
        self.edge_att = edge_att
        Event.__init__(self)

    def _calculate_state_variable(self, patch, edges):
        return sum([e[self.edge_att] for e in edges])

    def perform(self, network):
        pass


class NetworkTestCase(unittest.TestCase):

    def setUp(self):
        self.compartments = ['a','b','c']
        self.patch_attributes = ['d','e','f']
        self.edge_attributes = ['g', 'h', 'i']
        self.network = Network(self.compartments, self.patch_attributes, self.edge_attributes)

    def test_prepare(self):
        self.network.add_nodes_from(range(1, 8))
        edges = [(1, 2), (1, 3), (1, 4), (2, 5), (2, 6), (3, 7)]
        self.network.add_edges_from(edges)

        self.network.prepare()

        for n in range(1, 8):
            self.assertItemsEqual(self.network.node[n].keys(), [Network.COMPARTMENTS, Network.ATTRIBUTES,
                                                                Network.EVENTS, Network.TOTAL_RATE])
            self.assertItemsEqual(self.network.node[n][Network.COMPARTMENTS].keys(), self.compartments)
            for c in self.compartments:
                self.assertEqual(self.network.node[n][Network.COMPARTMENTS][c], 0)
            self.assertItemsEqual(self.network.node[n][Network.ATTRIBUTES].keys(), self.patch_attributes)
            for a in self.patch_attributes:
                self.assertEqual(self.network.node[n][Network.ATTRIBUTES][a], 0.0)

        for (u, v) in edges:
            data = self.network.edge[u][v]
            self.assertItemsEqual(data.keys(), self.edge_attributes)
            for a in self.edge_attributes:
                self.assertEqual(data[a], 0.0)

    def test_get_compartment_value(self):
        self.network.add_node(1)
        self.network.prepare()
        self.network.node[1][Network.COMPARTMENTS][self.compartments[0]] = 99
        self.assertEqual(self.network.get_compartment_value(1, self.compartments[0]), 99)

    def test_get_attribute_value(self):
        self.network.add_node(1)
        self.network.prepare()
        self.network.node[1][Network.ATTRIBUTES][self.patch_attributes[0]] = 99
        self.assertEqual(self.network.get_attribute_value(1, self.patch_attributes[0]), 99)

    def test_attach_event(self):
        self.network.add_node(1)
        self.network.prepare()

        e = NAEvent1(self.compartments[0], self.compartments[1])

        self.network.attach_event(e, 1)
        self.assertEqual(len(self.network.node[1][Network.EVENTS]), 1)
        self.assertItemsEqual(self.network.node[1][Network.EVENTS], [e])

    def test_update_rates_at_patch(self):
        self.network.add_node(1)
        self.network.prepare()

        e1 = NAEvent1(self.compartments[0], self.compartments[1])
        e2 = NAEvent1(self.compartments[1], self.compartments[2])

        self.network.attach_event(e1, 1)
        self.network.attach_event(e2, 1)

        self.network.node[1][Network.COMPARTMENTS][self.compartments[0]] = 1
        self.network.update_rates_at_patch(1)
        self.assertEqual(self.network.node[1][Network.TOTAL_RATE], e1.rate() + e2.rate())
        self.assertEqual(e1.rate(), e1.REACTION_PARAMETER * 1)
        self.assertEqual(e2.rate(), e1.REACTION_PARAMETER * 0)
        self.network.node[1][Network.COMPARTMENTS][self.compartments[1]] = 2
        self.network.update_rates_at_patch(1)
        self.assertEqual(self.network.node[1][Network.TOTAL_RATE], e1.rate() + e2.rate())
        self.assertEqual(e1.rate(), e1.REACTION_PARAMETER * 1)
        self.assertEqual(e2.rate(), e1.REACTION_PARAMETER * 2)

    def test_update_patch(self):
        self.network.add_nodes_from([1, 2, 3])
        self.network.prepare()

        NAEvent1.REACTION_PARAMETER = 0.1
        NAEvent2.REACTION_PARAMETER = 0.3

        for i in [1, 2, 3]:
            e1 = NAEvent1(self.compartments[0], self.compartments[1])
            e2 = NAEvent2(self.compartments[1], self.compartments[2])
            self.network.attach_event(e1, i)
            self.network.attach_event(e2, i)

        self.network.update_patch(1, {self.compartments[0]: 1, self.compartments[1]: 2},
                                  {self.patch_attributes[0]: 5, self.patch_attributes[1]: 6})
        self.network.update_patch(2, {self.compartments[0]: 3, self.compartments[1]: 4},
                                  {self.patch_attributes[0]: 7, self.patch_attributes[1]: 8})

        self.assertEqual(self.network.get_compartment_value(1, self.compartments[0]), 1)
        self.assertEqual(self.network.get_compartment_value(1, self.compartments[1]), 2)
        self.assertEqual(self.network.get_compartment_value(1, self.compartments[2]), 0)
        self.assertEqual(self.network.get_compartment_value(2, self.compartments[0]), 3)
        self.assertEqual(self.network.get_compartment_value(2, self.compartments[1]), 4)
        self.assertEqual(self.network.get_compartment_value(2, self.compartments[2]), 0)

        self.assertEqual(self.network.get_attribute_value(1, self.patch_attributes[0]), 5)
        self.assertEqual(self.network.get_attribute_value(1, self.patch_attributes[1]), 6)
        self.assertEqual(self.network.get_attribute_value(1, self.patch_attributes[2]), 0)
        self.assertEqual(self.network.get_attribute_value(2, self.patch_attributes[0]), 7)
        self.assertEqual(self.network.get_attribute_value(2, self.patch_attributes[1]), 8)
        self.assertEqual(self.network.get_attribute_value(2, self.patch_attributes[2]), 0)

        self.assertEqual(self.network.node[1][Network.TOTAL_RATE],
                         (NAEvent1.REACTION_PARAMETER * 1) + (NAEvent2.REACTION_PARAMETER * 1.5 * 2))
        self.assertEqual(self.network.node[2][Network.TOTAL_RATE],
                         (NAEvent1.REACTION_PARAMETER * 3) + (NAEvent2.REACTION_PARAMETER * 1.5 * 4))

    def test_update_edge(self):
        self.network.add_nodes_from([1, 2, 3])
        self.network.add_edges_from([(1,2), (2,3)])
        self.network.prepare()

        NAEvent3.REACTION_PARAMETER = 0.1

        for i in [1, 2, 3]:
            e1 = NAEvent3(self.edge_attributes[0])
            self.network.attach_event(e1, i)

        self.network.update_edge(1,2,{self.edge_attributes[0]: 2.3})

        self.assertEqual(self.network.node[1][Network.TOTAL_RATE],
                          (NAEvent3.REACTION_PARAMETER * 2.3))
        self.assertEqual(self.network.node[2][Network.TOTAL_RATE],
                         (NAEvent3.REACTION_PARAMETER * 2.3))


class TypedNetworkTestCase(unittest.TestCase):

    def setUp(self):
        self.patch_types = ['alpha', 'beta', 'gamma']
        self.compartments = ['a', 'b', 'c']
        self.patch_attributes = ['d', 'e', 'f']
        self.edge_attributes = ['g', 'h', 'i']
        self.network = TypedNetwork(self.compartments, self.patch_attributes, self.edge_attributes)

    def test_set_get_patch_type(self):
        self.network.add_nodes_from([1,2,3])
        self.network.set_patch_type(1, self.patch_types[0])
        self.network.set_patch_type(2, self.patch_types[0])
        self.network.set_patch_type(3, self.patch_types[1])

        self.network.prepare()
        self.assertItemsEqual(self.network.get_patches_by_type(self.patch_types[0]), [1, 2])
        self.assertItemsEqual(self.network.get_patches_by_type(self.patch_types[1]), [3])
        self.assertFalse(self.network.get_patches_by_type(self.patch_types[2]))


if __name__ == '__main__':
    unittest.main()
