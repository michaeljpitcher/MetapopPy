import unittest
from tbmetapoppy.events.cell_death import *
from metapoppy.network import Network


class CellDeathTestCase(unittest.TestCase):

    def setUp(self):
        self.event = CellDeath(TBDynamics.MACROPHAGE_RESTING)
        self.event_mi = CellDeath(TBDynamics.MACROPHAGE_INFECTED)
        self.event_dcm = CellDeath(TBDynamics.DENDRITIC_CELL_MATURE)

        self.rp = 0.1

        for e in [self.event, self.event_mi, self.event_dcm]:
            e.set_reaction_parameter(self.rp)

        self.network = Network(TBDynamics.TB_COMPARTMENTS, [], [])
        self.network.add_node(1)
        self.network.prepare()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {TBDynamics.MACROPHAGE_RESTING: 1})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.rp * 1)
        self.network.update_patch(1, {TBDynamics.MACROPHAGE_RESTING: 3})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.rp * 4)

    def test_perform(self):
        self.network.update_patch(1, {TBDynamics.MACROPHAGE_RESTING: 1})
        self.event.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBDynamics.MACROPHAGE_RESTING), 0)

        # Even split - 10/2
        self.network.update_patch(1, {TBDynamics.MACROPHAGE_INFECTED: 1,
                                      TBDynamics.BACTERIUM_INTRACELLULAR_MACROPHAGE: 10})
        self.event_mi.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBDynamics.MACROPHAGE_INFECTED), 0)
        self.assertEqual(self.network.get_compartment_value(1, TBDynamics.BACTERIUM_INTRACELLULAR_MACROPHAGE), 0)
        self.assertEqual(self.network.get_compartment_value(1, TBDynamics.BACTERIUM_EXTRACELLULAR_DORMANT), 10)

        # Round down - 10/3
        self.network.update_patch(1, {TBDynamics.BACTERIUM_EXTRACELLULAR_DORMANT: -10,
                                      TBDynamics.MACROPHAGE_INFECTED: 3,
                                      TBDynamics.BACTERIUM_INTRACELLULAR_MACROPHAGE: 10})
        self.event_mi.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBDynamics.MACROPHAGE_INFECTED), 2)
        self.assertEqual(self.network.get_compartment_value(1, TBDynamics.BACTERIUM_INTRACELLULAR_MACROPHAGE), 7)
        self.assertEqual(self.network.get_compartment_value(1, TBDynamics.BACTERIUM_EXTRACELLULAR_DORMANT), 3)

        # Round up - 10/4
        self.network.update_patch(1, {TBDynamics.BACTERIUM_EXTRACELLULAR_DORMANT: -3,
                                      TBDynamics.MACROPHAGE_INFECTED: 2,
                                      TBDynamics.BACTERIUM_INTRACELLULAR_MACROPHAGE: 3})
        self.event_mi.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBDynamics.MACROPHAGE_INFECTED), 3)
        self.assertEqual(self.network.get_compartment_value(1, TBDynamics.BACTERIUM_INTRACELLULAR_MACROPHAGE), 7)
        self.assertEqual(self.network.get_compartment_value(1, TBDynamics.BACTERIUM_EXTRACELLULAR_DORMANT), 3)

        self.network.update_patch(1, {TBDynamics.BACTERIUM_EXTRACELLULAR_DORMANT: -3,
                                      TBDynamics.MACROPHAGE_INFECTED: -3,
                                      TBDynamics.BACTERIUM_INTRACELLULAR_MACROPHAGE: -7,
                                      TBDynamics.DENDRITIC_CELL_MATURE: 3,
                                      TBDynamics.BACTERIUM_INTRACELLULAR_DENDRITIC: 3})
        self.event_dcm.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBDynamics.DENDRITIC_CELL_MATURE), 2)
        self.assertEqual(self.network.get_compartment_value(1, TBDynamics.BACTERIUM_INTRACELLULAR_DENDRITIC), 2)


class MacrophageBurstingTestCase(unittest.TestCase):

    def setUp(self):
        self.event = MacrophageBursting()
        self.rp = 0.1
        self.event.set_reaction_parameter(self.rp)
        self.network = Network(TBDynamics.TB_COMPARTMENTS, [], [])
        self.network.add_node(1)
        self.network.prepare()

    def test_rate(self):
        self.event.set_parameters(2, 20)
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {TBDynamics.MACROPHAGE_INFECTED: 3,
                                      TBDynamics.BACTERIUM_INTRACELLULAR_MACROPHAGE: 10})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1),
                         self.rp * 3.0 * (10.0 ** 2) / (10**2 + (20 * 3) ** 2))


class TCellDestroysMacrophageTestCase(unittest.TestCase):

    def setUp(self):
        self.event = TCellDestroysMacrophage()
        self.rp = 0.1
        self.event.set_reaction_parameter(self.rp)
        self.network = Network(TBDynamics.TB_COMPARTMENTS, [], [])
        self.network.add_node(1)
        self.network.prepare()

    def test_rate(self):
        self.event.set_parameters(10)
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {TBDynamics.MACROPHAGE_INFECTED: 3, TBDynamics.T_CELL_ACTIVATED: 7})
        self.assertAlmostEqual(self.event.calculate_rate_at_patch(self.network, 1), self.rp * 3.0 * 7.0 / (7 + 10))


class MacrophageDestroysBacteriumTestCase(unittest.TestCase):

    def setUp(self):
        self.event = MacrophageDestroysBacterium(TBDynamics.MACROPHAGE_RESTING,
                                                 TBDynamics.BACTERIUM_EXTRACELLULAR_REPLICATING)
        self.rp = 0.1
        self.event.set_reaction_parameter(self.rp)
        self.network = Network(TBDynamics.TB_COMPARTMENTS, [], [])
        self.network.add_node(1)
        self.network.prepare()

    def test_rate(self):
        self.event.set_parameters(100)
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {TBDynamics.MACROPHAGE_RESTING: 3,
                                      TBDynamics.BACTERIUM_EXTRACELLULAR_REPLICATING: 7})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.rp * 7.0 * 3.0 / (3 + 100))


if __name__ == '__main__':
    unittest.main()
