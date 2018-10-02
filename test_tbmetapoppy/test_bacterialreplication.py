import unittest
from tbmetapoppy import *
from metapoppy.network import Network


class ExtracellularBacterialReplicationTestCase(unittest.TestCase):

    def setUp(self):
        self.event = Replication(BACTERIUM_EXTRACELLULAR_REPLICATING)
        self.rp = 0.1
        self.event.set_reaction_parameter(self.rp)
        self.network = Network(TB_COMPARTMENTS, [], [])
        self.network.add_node(1)
        self.network.prepare()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {BACTERIUM_EXTRACELLULAR_REPLICATING: 1})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.rp * 1)
        self.network.update_patch(1, {BACTERIUM_EXTRACELLULAR_REPLICATING: 3})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.rp * 4)

    def test_perform(self):
        self.network.update_patch(1, {BACTERIUM_EXTRACELLULAR_REPLICATING: 1})
        self.event.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_EXTRACELLULAR_REPLICATING), 2)


class IntracellularBacterialReplicationTestCase(unittest.TestCase):

    def setUp(self):
        self.event = IntracellularBacterialReplication()
        self.rp = 0.1
        self.sig = 2
        self.hs = 30
        self.event.set_parameters(self.rp, self.sig, self.hs)
        self.network = Network(TB_COMPARTMENTS, [], [])
        self.network.add_node(1)
        self.network.prepare()

    def test_rate(self):

        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {BACTERIUM_INTRACELLULAR_MACROPHAGE: 18,
                                      MACROPHAGE_INFECTED: 5})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1),
                         0.1 * 18 * (1 - (1.0 * (18 ** self.sig) / (18 ** self.sig + (5 * self.hs) ** self.sig))))

    def test_perform(self):
        self.network.update_patch(1, {MACROPHAGE_INFECTED: 1,
                                      BACTERIUM_INTRACELLULAR_MACROPHAGE: 10})
        self.event.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_INTRACELLULAR_MACROPHAGE), 11)


if __name__ == '__main__':
    unittest.main()
