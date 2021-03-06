import unittest
from metapoppy import *
import numpy


class NetworkTestCase(unittest.TestCase):

    def setUp(self):
        self.compartments = ['a','b','c']
        self.patch_attributes = ['d','e','f']
        self.edge_attributes = ['g', 'h', 'i']
        self.network = Environment(self.compartments, self.patch_attributes, self.edge_attributes)

    def test_reset(self):
        self.network.add_nodes_from(range(1, 8))
        edges = [(1, 2), (1, 3), (1, 4), (2, 5), (2, 6), (3, 7)]
        self.network.add_edges_from(edges)

        self.network.reset()

        for n in range(1, 8):
            self.assertItemsEqual(self.network.node[n].keys(), [Environment.COMPARTMENTS, Environment.ATTRIBUTES])
            self.assertItemsEqual(self.network.node[n][Environment.COMPARTMENTS].keys(), self.compartments)
            for c in self.compartments:
                self.assertEqual(self.network.node[n][Environment.COMPARTMENTS][c], 0)
            self.assertItemsEqual(self.network.node[n][Environment.ATTRIBUTES].keys(), self.patch_attributes)
            for a in self.patch_attributes:
                self.assertEqual(self.network.node[n][Environment.ATTRIBUTES][a], 0.0)

        for (u, v) in edges:
            data = self.network.get_edge_data(u, v)
            self.assertItemsEqual(data.keys(), self.edge_attributes)
            for a in self.edge_attributes:
                self.assertEqual(data[a], 0.0)

    def test_get_compartment_value(self):
        self.network.add_node(1)
        self.network.reset()
        self.network.node[1][Environment.COMPARTMENTS][self.compartments[0]] = 99
        self.network.node[1][Environment.COMPARTMENTS][self.compartments[1]] = 2
        self.assertEqual(self.network.get_compartment_value(1, self.compartments[0]), 99)

        # Multiple compartments
        self.assertEqual(self.network.get_compartment_value(1, self.compartments[0:2]), 101)

    def test_get_attribute_value(self):
        self.network.add_node(1)
        self.network.reset()
        self.network.node[1][Environment.ATTRIBUTES][self.patch_attributes[0]] = 99
        self.network.node[1][Environment.ATTRIBUTES][self.patch_attributes[1]] = 2
        self.assertEqual(self.network.get_attribute_value(1, self.patch_attributes[0]), 99)

        # Multiple attributes
        self.assertEqual(self.network.get_attribute_value(1, self.patch_attributes[0:2]), 101)

    def test_update_patch(self):
        # With handler
        self.check_value = []

        # When an update is performed store values used in a variable
        def patch_handler(a,b,c):
            self.check_value.append([a,b,c])


        self.network.set_handlers(lambda a, b, c: patch_handler(a, b, c), None)

        self.network.add_node(1)
        self.network.reset()

        self.network.update_patch(1, {self.compartments[0]: 1, self.compartments[1]: 2},
                                  {self.patch_attributes[0]: 3, self.patch_attributes[1]: 4})

        # Check the values have been updated
        self.assertEqual(self.network.get_compartment_value(1, self.compartments[0]), 1)
        self.assertEqual(self.network.get_compartment_value(1, self.compartments[1]), 2)
        self.assertEqual(self.network.get_attribute_value(1, self.patch_attributes[0]), 3)
        self.assertEqual(self.network.get_attribute_value(1, self.patch_attributes[1]), 4)

        # Check the values were passed into the handler
        self.assertEqual(self.check_value[0][0], 1)
        self.assertItemsEqual(self.check_value[0][1], [self.compartments[0], self.compartments[1]])
        self.assertItemsEqual(self.check_value[0][2], [self.patch_attributes[0], self.patch_attributes[1]])

    def test_update_edge(self):
        self.check_value = []

        # When an update is performed store values used in a variable
        def edge_handler(a, b, c):
            self.check_value.append([a, b, c])

        self.network.add_node(1)
        self.network.add_node(2)
        self.network.add_edge(1, 2)
        self.network.reset()
        self.network.set_handlers(None, lambda a, b, c: edge_handler(a, b, c))

        self.network.update_edge(1, 2, {self.edge_attributes[0]: 1, self.edge_attributes[1]: 2})

        self.assertEqual(self.check_value[0][0], 1)
        self.assertEqual(self.check_value[0][1], 2)
        self.assertItemsEqual(self.check_value[0][2], [self.edge_attributes[0], self.edge_attributes[1]])


class TypedNetworkTestCase(unittest.TestCase):

    def setUp(self):
        self.patch_types = ['alpha', 'beta', 'gamma']
        self.compartments = ['a', 'b', 'c']
        self.patch_attributes = {self.patch_types[0]: ['d', 'e'], self.patch_types[1]: ['f']}
        self.edge_attributes = ['g', 'h', 'i']
        self.network = TypedEnvironment(self.compartments, self.patch_attributes, self.edge_attributes)

    def test_set_get_patch_type(self):
        self.network.add_nodes_from([1,2,3])
        self.network.set_patch_type(1, self.patch_types[0])
        self.network.set_patch_type(2, self.patch_types[0])
        self.network.set_patch_type(3, self.patch_types[1])

        self.network.reset()
        self.assertItemsEqual(self.network.get_patches_by_type(self.patch_types[0]), [1, 2])
        self.assertItemsEqual(self.network.get_patches_by_type(self.patch_types[1]), [3])
        self.assertFalse(self.network.get_patches_by_type(self.patch_types[2]))

        self.assertItemsEqual(self.network.get_patches_by_type(self.patch_types[0], data=True),
                              [(d, self.network.node[d]) for d in [1, 2]])
        self.assertItemsEqual(self.network.get_patches_by_type(self.patch_types[1], data=True),
                              [(d, self.network.node[d]) for d in [3]])

        # Correct attributes
        self.network.add_node(4)
        self.network.set_patch_type(4, self.patch_types[2])
        self.network.reset()
        for _,d in self.network.nodes(data=True):
            if d[TypedEnvironment.PATCH_TYPE] == self.patch_types[0]:
                self.assertItemsEqual(d[Environment.ATTRIBUTES], self.patch_attributes[self.patch_types[0]])
            elif d[TypedEnvironment.PATCH_TYPE] == self.patch_types[1]:
                self.assertItemsEqual(d[Environment.ATTRIBUTES], self.patch_attributes[self.patch_types[1]])
            elif d[TypedEnvironment.PATCH_TYPE] == self.patch_types[2]:
                self.assertFalse(d[Environment.ATTRIBUTES])


if __name__ == '__main__':
    unittest.main()
