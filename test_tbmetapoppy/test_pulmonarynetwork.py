import unittest
from tbmetapoppy import *


class PulmonaryNetworkTestCase(unittest.TestCase):

    def setUp(self):
        self.boundary = [(0,5), (0,10), (10,10), (10,0), (0,0)]
        self.network_config = {PulmonaryNetwork.TOPOLOGY: PulmonaryNetwork.SPACE_FILLING_TREE_2D,
                               PulmonaryNetwork.BOUNDARY: self.boundary,
                               PulmonaryNetwork.LENGTH_DIVISOR: 2,
                               PulmonaryNetwork.MINIMUM_AREA: 1}
        self.compartments = ['a', 'b', 'c']
        self.tree_network = PulmonaryNetwork(self.network_config, self.compartments)
        self.tree_network.prepare()

    def test_2d_space_filling_tree_initialise(self):
        print self.tree_network.get_patches_by_type(PulmonaryNetwork.ALVEOLAR_PATCH)
        print self.tree_network._alveolar_positions
        print self.tree_network.edges(data=True)




if __name__ == '__main__':
    unittest.main()
