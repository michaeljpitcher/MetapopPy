import unittest
from tbmetapoppy import *


class PulmonaryNetworkTestCase(unittest.TestCase):

    def setUp(self):
        self.boundary = [(0,5), (0,10), (10,10), (10,0), (0,0)]
        self.network_config = {PulmonaryNetwork.TOPOLOGY: PulmonaryNetwork.SPACE_FILLING_TREE_2D,
                               PulmonaryNetwork.BOUNDARY: self.boundary,
                               PulmonaryNetwork.LENGTH_DIVISOR: 2,
                               PulmonaryNetwork.MINIMUM_AREA: 6}
        self.compartments = ['a', 'b', 'c']
        self.tree_network = PulmonaryNetwork(self.network_config, self.compartments)
        self.tree_network.prepare()

    def test_2d_space_filling_tree_initialise(self):
        self.assertEqual(len(self.tree_network.nodes()), 17)
        self.assertEqual(len(self.tree_network.get_patches_by_type(PulmonaryNetwork.ALVEOLAR_PATCH)), 16)
        self.assertEqual(len(self.tree_network.get_patches_by_type(PulmonaryNetwork.LYMPH_PATCH)), 1)

        positions = [1.25, 3.75, 6.25, 8.75]
        e = []
        for q in positions:
            for w in positions:
                e.append((q,w))
        self.assertItemsEqual(e, self.tree_network._alveolar_positions.values())
        self.assertEqual(len(self.tree_network.edges()), 16)
        self.assertEqual(len(self.tree_network.edges([PulmonaryNetwork.LYMPH_PATCH])), 16)

    def test_pulmonary_attribute_seeding(self):
        seeding = self.tree_network.pulmonary_attribute_seeding(3,2,1)

        self.assertEqual(sum(a[PulmonaryNetwork.VENTILATION] for a in seeding.values()), 1)
        self.assertEqual(sum(a[PulmonaryNetwork.PERFUSION] for a in seeding.values()), 1)
        self.assertEqual(sum(a[PulmonaryNetwork.DRAINAGE] for a in seeding.values()), 1)

        # Check all drainage values are the same (since we've set the skew to 1)
        self.assertEqual(len(set([a[PulmonaryNetwork.DRAINAGE] for a in seeding.values()])), 1)

        for d in seeding.values():
            self.assertEqual(d[PulmonaryNetwork.OXYGEN_TENSION],
                             d[PulmonaryNetwork.VENTILATION] / d[PulmonaryNetwork.PERFUSION])


if __name__ == '__main__':
    unittest.main()
