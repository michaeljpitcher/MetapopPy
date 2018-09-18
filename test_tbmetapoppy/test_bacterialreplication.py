import unittest
from tbmetapoppy import *


class ExtracellularBacterialReplicationTestCase(unittest.TestCase):

    def setUp(self):
        self.event = ExtracellularBacterialReplication(TBDynamics.BACTERIUM_EXTRACELLULAR_REPLICATING)
        self.rp = 0.1
        self.event.set_reaction_parameter(self.rp)
        self.network = Network(TBDynamics.TB_COMPARTMENTS, [], [])
        self.network.add_node(1)
        self.network.prepare()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {TBDynamics.BACTERIUM_EXTRACELLULAR_REPLICATING: 1})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.rp * 1)
        self.network.update_patch(1, {TBDynamics.BACTERIUM_EXTRACELLULAR_REPLICATING: 3})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.rp * 4)

    def test_perform(self):
        self.network.update_patch(1, {TBDynamics.BACTERIUM_EXTRACELLULAR_REPLICATING: 1})
        self.event.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBDynamics.BACTERIUM_EXTRACELLULAR_REPLICATING), 2)


if __name__ == '__main__':
    unittest.main()
