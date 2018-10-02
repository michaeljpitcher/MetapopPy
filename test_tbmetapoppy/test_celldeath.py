import unittest
from tbmetapoppy import *
from metapoppy.network import Network


class CellDeathTestCase(unittest.TestCase):

    def setUp(self):
        self.event = CellDeath(MACROPHAGE_RESTING)
        self.event_mi = CellDeath(MACROPHAGE_INFECTED)
        self.event_dcm = CellDeath(DENDRITIC_CELL_MATURE)

        self.rp = 0.1

        for e in [self.event, self.event_mi, self.event_dcm]:
            e.set_reaction_parameter(self.rp)

        self.network = Network(TB_COMPARTMENTS, [], [])
        self.network.add_node(1)
        self.network.prepare()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {MACROPHAGE_RESTING: 1})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.rp * 1)
        self.network.update_patch(1, {MACROPHAGE_RESTING: 3})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.rp * 4)

    def test_perform(self):
        self.network.update_patch(1, {MACROPHAGE_RESTING: 1})
        self.event.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, MACROPHAGE_RESTING), 0)

        # Even split - 10/2
        self.network.update_patch(1, {MACROPHAGE_INFECTED: 1,
                                      BACTERIUM_INTRACELLULAR_MACROPHAGE: 10})
        self.event_mi.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, MACROPHAGE_INFECTED), 0)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_INTRACELLULAR_MACROPHAGE), 0)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_EXTRACELLULAR_DORMANT), 10)

        # Round down - 10/3
        self.network.update_patch(1, {BACTERIUM_EXTRACELLULAR_DORMANT: -10,
                                      MACROPHAGE_INFECTED: 3,
                                      BACTERIUM_INTRACELLULAR_MACROPHAGE: 10})
        self.event_mi.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, MACROPHAGE_INFECTED), 2)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_INTRACELLULAR_MACROPHAGE), 7)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_EXTRACELLULAR_DORMANT), 3)

        # Round up - 10/4
        self.network.update_patch(1, {BACTERIUM_EXTRACELLULAR_DORMANT: -3,
                                      MACROPHAGE_INFECTED: 2,
                                      BACTERIUM_INTRACELLULAR_MACROPHAGE: 3})
        self.event_mi.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, MACROPHAGE_INFECTED), 3)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_INTRACELLULAR_MACROPHAGE), 7)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_EXTRACELLULAR_DORMANT), 3)

        self.network.update_patch(1, {BACTERIUM_EXTRACELLULAR_DORMANT: -3,
                                      MACROPHAGE_INFECTED: -3,
                                      BACTERIUM_INTRACELLULAR_MACROPHAGE: -7,
                                      DENDRITIC_CELL_MATURE: 3,
                                      BACTERIUM_INTRACELLULAR_DENDRITIC: 3})
        self.event_dcm.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, DENDRITIC_CELL_MATURE), 2)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_INTRACELLULAR_DENDRITIC), 2)


class MacrophageBurstingTestCase(unittest.TestCase):

    def setUp(self):
        self.event = MacrophageBursting()
        self.rp = 0.1
        self.sig = 2
        self.cap = 20
        self.event.set_parameters(self.rp, self.sig, self.cap)
        self.network = Network(TB_COMPARTMENTS, [], [])
        self.network.add_node(1)
        self.network.prepare()

    def test_rate(self):

        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {MACROPHAGE_INFECTED: 3,
                                      BACTERIUM_INTRACELLULAR_MACROPHAGE: 10})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1),
                         self.rp * 3.0 * (10.0 ** self.sig) / (10**self.sig + (self.cap * 3) ** self.sig))


class TCellDestroysMacrophageTestCase(unittest.TestCase):

    def setUp(self):
        self.event = TCellDestroysMacrophage()
        self.rp = 0.1
        self.hs = 10
        self.event.set_parameters(self.rp, self.hs)
        self.network = Network(TB_COMPARTMENTS, [], [])
        self.network.add_node(1)
        self.network.prepare()

    def test_rate(self):

        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {MACROPHAGE_INFECTED: 3, T_CELL_ACTIVATED: 7})
        self.assertAlmostEqual(self.event.calculate_rate_at_patch(self.network, 1), self.rp * 3.0 * 7.0 / (7 + self.hs))


class MacrophageDestroysBacteriumTestCase(unittest.TestCase):

    def setUp(self):
        self.event = MacrophageDestroysBacterium(MACROPHAGE_RESTING,
                                                 BACTERIUM_EXTRACELLULAR_REPLICATING)
        self.rp = 0.1
        self.hs = 100
        self.event.set_parameters(self.rp, self.hs)
        self.network = Network(TB_COMPARTMENTS, [], [])
        self.network.add_node(1)
        self.network.prepare()

    def test_rate(self):

        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {MACROPHAGE_RESTING: 3,
                                      BACTERIUM_EXTRACELLULAR_REPLICATING: 7})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.rp * 7.0 * 3.0 / (3 + self.hs))


if __name__ == '__main__':
    unittest.main()
