import unittest
from tbmetapoppy import *
from metapoppy.network import Network

TEST_RATE_DEATH_MR = 'test_rate_death_mr'
TEST_RATE_DEATH_MI = 'test_rate_death_mi'
TEST_RATE_DEATH_DCM = 'test_rate_death_dcm'
TEST_RATE_BURST = 'test_rate_burst'
TEST_SIGMOID = 'test_sigmoid'
TEST_CAPACITY = 'test_capacity'
TEST_RATE_TA_KILL = 'test_rate_ta_kill'
TEST_HALFSAT_TA_KILL = 'test_halfsat_ta_kill'
TEST_RATE_MACBAC_KILL = 'test_rate_macbac_kill'
TEST_HALFSAT_MACBAC_KILL = 'test_halfsat_macbac_kill'


class CellDeathTestCase(unittest.TestCase):

    def setUp(self):
        self.event = CellDeath(TEST_RATE_DEATH_MR, MACROPHAGE_RESTING)
        self.event_mi = CellDeath(TEST_RATE_DEATH_MI, MACROPHAGE_INFECTED)
        self.event_dcm = CellDeath(TEST_RATE_DEATH_DCM, DENDRITIC_CELL_MATURE)

        self.params = {TEST_RATE_DEATH_MR: 0.1, TEST_RATE_DEATH_MI: 0.2, TEST_RATE_DEATH_DCM: 0.3}
        for e in [self.event, self.event_mi, self.event_dcm]:
            e.set_parameters(self.params)

        self.network = Network(TB_COMPARTMENTS, [], [])
        self.network.add_node(1)
        self.network.prepare()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {MACROPHAGE_RESTING: 1, MACROPHAGE_INFECTED: 2})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.params[TEST_RATE_DEATH_MR] * 1)
        self.assertEqual(self.event_mi.calculate_rate_at_patch(self.network, 1), self.params[TEST_RATE_DEATH_MI] * 2)
        self.network.update_patch(1, {MACROPHAGE_RESTING: 3})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.params[TEST_RATE_DEATH_MR] * 4)

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
        self.assertEqual(self.network.get_compartment_value(1, CASEUM), 1)

        # Round down - 10/3
        self.network.update_patch(1, {BACTERIUM_EXTRACELLULAR_DORMANT: -10,
                                      MACROPHAGE_INFECTED: 3,
                                      BACTERIUM_INTRACELLULAR_MACROPHAGE: 10})
        self.event_mi.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, MACROPHAGE_INFECTED), 2)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_INTRACELLULAR_MACROPHAGE), 7)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_EXTRACELLULAR_DORMANT), 3)
        self.assertEqual(self.network.get_compartment_value(1, CASEUM), 2)

        # Round up - 10/4
        self.network.update_patch(1, {BACTERIUM_EXTRACELLULAR_DORMANT: -3,
                                      MACROPHAGE_INFECTED: 2,
                                      BACTERIUM_INTRACELLULAR_MACROPHAGE: 3})
        self.event_mi.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, MACROPHAGE_INFECTED), 3)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_INTRACELLULAR_MACROPHAGE), 7)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_EXTRACELLULAR_DORMANT), 3)
        self.assertEqual(self.network.get_compartment_value(1, CASEUM), 3)

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
        self.event = MacrophageBursting(TEST_RATE_BURST, TEST_SIGMOID, TEST_CAPACITY)
        self.params = {TEST_RATE_BURST: 0.1, TEST_SIGMOID: 2, TEST_CAPACITY: 20}
        self.event.set_parameters(self.params)
        self.network = Network(TB_COMPARTMENTS, [], [])
        self.network.add_node(1)
        self.network.prepare()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {MACROPHAGE_INFECTED: 3,
                                      BACTERIUM_INTRACELLULAR_MACROPHAGE: 10})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1),
                         self.params[TEST_RATE_BURST] * 3.0 * (10.0 ** self.params[TEST_SIGMOID]) /
                         (10**self.params[TEST_SIGMOID] +
                          (self.params[TEST_CAPACITY] * 3) ** self.params[TEST_SIGMOID]))


class TCellDestroysMacrophageTestCase(unittest.TestCase):

    def setUp(self):
        self.event = TCellDestroysMacrophage(TEST_RATE_TA_KILL, TEST_HALFSAT_TA_KILL)
        self.params = {TEST_RATE_TA_KILL: 0.1, TEST_HALFSAT_TA_KILL: 10}
        self.event.set_parameters(self.params)
        self.network = Network(TB_COMPARTMENTS, [], [])
        self.network.add_node(1)
        self.network.prepare()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {MACROPHAGE_INFECTED: 3, T_CELL_ACTIVATED: 7})
        self.assertAlmostEqual(self.event.calculate_rate_at_patch(self.network, 1),
                               self.params[TEST_RATE_TA_KILL] * 3.0 * 7.0 / (7 + self.params[TEST_HALFSAT_TA_KILL]))


class MacrophageDestroysBacteriumTestCase(unittest.TestCase):

    def setUp(self):
        self.event = MacrophageDestroysBacterium(TEST_RATE_MACBAC_KILL, MACROPHAGE_RESTING,
                                                 TEST_HALFSAT_MACBAC_KILL)
        self.params = {TEST_RATE_MACBAC_KILL: 0.1, TEST_HALFSAT_MACBAC_KILL: 100}
        self.event.set_parameters(self.params)
        self.network = Network(TB_COMPARTMENTS, [], [])
        self.network.add_node(1)
        self.network.prepare()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {MACROPHAGE_RESTING: 3,
                                      BACTERIUM_EXTRACELLULAR_REPLICATING: 7})
        self.assertAlmostEqual(self.event.calculate_rate_at_patch(self.network, 1),
                         self.params[TEST_RATE_MACBAC_KILL] * 3.0 * 7.0 / (7 + self.params[TEST_HALFSAT_MACBAC_KILL]))
        self.network.update_patch(1, {BACTERIUM_EXTRACELLULAR_DORMANT: 11})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1),
                         self.params[TEST_RATE_MACBAC_KILL] * 3.0 * 18.0 / (18 + self.params[TEST_HALFSAT_MACBAC_KILL]))


if __name__ == '__main__':
    unittest.main()
