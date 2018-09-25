import unittest
from tbmetapoppy import *


class TBDynamicsTestCase(unittest.TestCase):

    def setUp(self):
        self.boundary = [(0, 5), (0, 10), (10, 10), (10, 0), (0, 0)]
        network_config = {PulmonaryNetwork.TOPOLOGY: PulmonaryNetwork.SPACE_FILLING_TREE_2D,
                          PulmonaryNetwork.BOUNDARY: self.boundary,
                          PulmonaryNetwork.LENGTH_DIVISOR: 2,
                          PulmonaryNetwork.MINIMUM_AREA: 6}
        self.dynamics = TBDynamics(network_config)

    def test_initialise(self):
        # Rows = number of events, cols = number of nodes
        self.assertEqual(self.dynamics._rate_table.shape, (32, 17))




if __name__ == '__main__':
    unittest.main()
