import unittest
from metapoppy import *
import numpy


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
            self.assertItemsEqual(self.network.node[n].keys(), [Network.COMPARTMENTS, Network.ATTRIBUTES])
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

    def test_update_patch(self):
        # With handler
        self.check_value = []

        # When an update is performed store values used in a variable
        def handler(a,b,c,d):
            self.check_value = [a,b,c,d]

        self.network.add_node(1)
        self.network.prepare(lambda a,b,c,d: handler(a,b,c,d))

        self.network.update_patch(1, {self.compartments[0]: 1, self.compartments[1]: 2},
                                  {self.patch_attributes[0]: 3, self.patch_attributes[1]: 4})

        # Check the values have been updated
        self.assertEqual(self.network.get_compartment_value(1, self.compartments[0]), 1)
        self.assertEqual(self.network.get_compartment_value(1, self.compartments[1]), 2)
        self.assertEqual(self.network.get_attribute_value(1, self.patch_attributes[0]), 3)
        self.assertEqual(self.network.get_attribute_value(1, self.patch_attributes[1]), 4)

        # Check the values were passed into the handler
        # print self.check_value
        self.assertEqual(self.check_value[0], 1)
        self.assertItemsEqual(self.check_value[1], [self.compartments[0], self.compartments[1]])
        self.assertItemsEqual(self.check_value[2], [self.patch_attributes[0], self.patch_attributes[1]])
        self.assertFalse(self.check_value[3])

    def test_update_edge(self):
        self.check_value = []

        # When an update is performed store values used in a variable
        def handler(a, b, c, d):
            self.check_value.append([a, b, c, d])

        self.network.add_node(1)
        self.network.add_node(2)
        self.network.add_edge(1, 2)
        self.network.prepare(lambda a, b, c, d: handler(a, b, c, d))

        self.network.update_edge(1, 2, {self.edge_attributes[0]: 1, self.edge_attributes[1]: 2})

        self.assertEqual(self.check_value[0][0], 1)
        self.assertFalse(self.check_value[0][1])
        self.assertFalse(self.check_value[0][2])
        self.assertItemsEqual(self.check_value[0][3], [self.edge_attributes[0], self.edge_attributes[1]])


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
