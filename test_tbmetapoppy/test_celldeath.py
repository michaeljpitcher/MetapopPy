import unittest
from tbmetapoppy import *
from metapoppy.network import MetapopulationNetwork



class CellDeathTestCase(unittest.TestCase):

    def setUp(self):
        self.event = CellDeath(TBPulmonaryNetwork.MACROPHAGE_RESTING)
        self.params = {'m_r_death_rate': 0.1}
        self.event.set_parameters(self.params)

        self.network = TBPulmonaryNetwork({TBPulmonaryNetwork.TOPOLOGY: None})
        self.network.add_node(1)
        self.network.set_patch_type(1, TBPulmonaryNetwork.LYMPH_PATCH)
        self.network.reset()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {TBPulmonaryNetwork.MACROPHAGE_RESTING: 1,
                                      TBPulmonaryNetwork.MACROPHAGE_INFECTED: 2})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.params['m_r_death_rate'] * 1)
        self.network.update_patch(1, {TBPulmonaryNetwork.MACROPHAGE_RESTING: 3})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.params['m_r_death_rate'] * 4)

    def test_perform(self):
        self.network.update_patch(1, {TBPulmonaryNetwork.MACROPHAGE_RESTING: 1})
        self.event.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.MACROPHAGE_RESTING), 0)

        self.network.update_patch(1, {TBPulmonaryNetwork.DENDRITIC_CELL_MATURE: 3,
                                      TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_DENDRITIC: 3})


class InfectedCellDeathTestCase(unittest.TestCase):

    def setUp(self):
        self.event_mi = InfectedCellDeath(TBPulmonaryNetwork.MACROPHAGE_INFECTED)
        self.event_dcm = InfectedCellDeath(TBPulmonaryNetwork.DENDRITIC_CELL_MATURE)

        self.params = {'m_i_death_rate': 0.1, 'm_i_death_percentage_bacteria_destroyed': 0.6, 'd_m_death_rate': 0.2,
                       'd_m_death_percentage_bacteria_destroyed': 0.0}
        self.event_mi.set_parameters(self.params)
        self.event_dcm.set_parameters(self.params)

        self.network = TBPulmonaryNetwork({TBPulmonaryNetwork.TOPOLOGY: None})
        self.network.add_node(1)
        self.network.set_patch_type(1, TBPulmonaryNetwork.LYMPH_PATCH)
        self.network.reset()

    def test_rate(self):
        self.assertFalse(self.event_mi.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {TBPulmonaryNetwork.MACROPHAGE_INFECTED: 2})
        self.assertEqual(self.event_mi.calculate_rate_at_patch(self.network, 1), self.params['m_i_death_rate'] * 2)
        self.network.update_patch(1, {TBPulmonaryNetwork.DENDRITIC_CELL_MATURE: 66})
        self.assertEqual(self.event_mi.calculate_rate_at_patch(self.network, 1), self.params['m_i_death_rate'] * 2)
        self.assertEqual(self.event_dcm.calculate_rate_at_patch(self.network, 1), self.params['d_m_death_rate'] * 66)

    def test_perform(self):
        self.network.update_patch(1, {TBPulmonaryNetwork.MACROPHAGE_INFECTED: 1,
                                      TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE: 10})
        self.event_mi.perform(self.network, 1)
        exp_kill = round(self.params['m_i_death_percentage_bacteria_destroyed'] * (10.0/1))
        exp_rel = round((1 - self.params['m_i_death_percentage_bacteria_destroyed']) * (10.0/1))
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.MACROPHAGE_INFECTED), 0)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE),
                         10 - exp_kill - exp_rel)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT), exp_rel)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.CASEUM), 1)

        self.network.reset()
        self.network.update_patch(1, {TBPulmonaryNetwork.MACROPHAGE_INFECTED: 4,
                                      TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE: 10})
        self.event_mi.perform(self.network, 1)
        exp_kill = round(self.params['m_i_death_percentage_bacteria_destroyed'] * (10.0 / 4))
        exp_rel = round((1 - self.params['m_i_death_percentage_bacteria_destroyed']) * (10.0 / 4))
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.MACROPHAGE_INFECTED), 3)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE),
                         10 - exp_kill - exp_rel)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT),
                         exp_rel)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.CASEUM), 1)

        self.network.reset()
        self.network.update_patch(1, {TBPulmonaryNetwork.DENDRITIC_CELL_MATURE: 4,
                                      TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_DENDRITIC: 4})
        self.event_dcm.perform(self.network, 1)
        exp_kill = round(self.params['d_m_death_percentage_bacteria_destroyed'] * (4.0 / 4))
        exp_rel = round((1 - self.params['d_m_death_percentage_bacteria_destroyed']) * (4.0 / 4))
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.DENDRITIC_CELL_MATURE), 3)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_DENDRITIC),
                         4 - exp_kill - exp_rel)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT),
                         exp_rel)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.CASEUM), 0)


class MacrophageBurstingTestCase(unittest.TestCase):

    def setUp(self):
        self.event = MacrophageBursting()
        self.params = {MacrophageBursting.BURSTING + RATE: 0.1, INTRACELLULAR_REPLICATION_SIGMOID: 2,
                       MACROPHAGE_CAPACITY: 33, 'macrophage_bursting_percentage_bacteria_destroyed': 0.0}
        self.event.set_parameters(self.params)
        self.network = TBPulmonaryNetwork({TBPulmonaryNetwork.TOPOLOGY: None})
        self.network.add_node(1)
        self.network.set_patch_type(1, TBPulmonaryNetwork.LYMPH_PATCH)
        self.network.reset()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {TBPulmonaryNetwork.MACROPHAGE_INFECTED: 3,
                                      TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE: 10})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1),
                         self.params[MacrophageBursting.BURSTING + RATE] * 3.0 * (10.0 ** self.params[INTRACELLULAR_REPLICATION_SIGMOID]) /
                         (10**self.params[INTRACELLULAR_REPLICATION_SIGMOID] +
                          (self.params[MACROPHAGE_CAPACITY] * 3) ** self.params[INTRACELLULAR_REPLICATION_SIGMOID]))


class TCellDestroysMacrophageTestCase(unittest.TestCase):

    def setUp(self):
        self.event = TCellDestroysMacrophage()
        self.params = {'t_cell_destroys_macrophage_rate': 0.1, 't_cell_destroys_macrophage_half_sat': 25, 't_cell_destroys_macrophage_percentage_bacteria_destroyed': 0.5}
        self.event.set_parameters(self.params)
        self.network = TBPulmonaryNetwork({TBPulmonaryNetwork.TOPOLOGY: None})
        self.network.add_node(1)
        self.network.set_patch_type(1, TBPulmonaryNetwork.LYMPH_PATCH)
        self.network.reset()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {TBPulmonaryNetwork.MACROPHAGE_INFECTED: 3, TBPulmonaryNetwork.T_CELL_ACTIVATED: 7})
        self.assertAlmostEqual(self.event.calculate_rate_at_patch(self.network, 1),
                               self.params['t_cell_destroys_macrophage_rate'] * 3.0 * 7.0 / (7 + self.params['t_cell_destroys_macrophage_half_sat']))


class MacrophageDestroysBacteriumTestCase(unittest.TestCase):

    def setUp(self):
        self.event = MacrophageDestroysBacterium(TBPulmonaryNetwork.MACROPHAGE_RESTING)
        self.params = {'m_r_destroys_bacterium_rate': 0.1, 'm_r_destroys_bacterium_half_sat': 100}
        self.event.set_parameters(self.params)
        self.network = TBPulmonaryNetwork({TBPulmonaryNetwork.TOPOLOGY: None})
        self.network.add_node(1)
        self.network.set_patch_type(1, TBPulmonaryNetwork.LYMPH_PATCH)
        self.network.reset()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {TBPulmonaryNetwork.MACROPHAGE_RESTING: 3,
                                      TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING: 7})
        self.assertAlmostEqual(self.event.calculate_rate_at_patch(self.network, 1),
                         self.params['m_r_destroys_bacterium_rate'] * 3.0 * 7.0 / (7 + self.params['m_r_destroys_bacterium_half_sat']))
        self.network.update_patch(1, {TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT: 11})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1),
                         self.params['m_r_destroys_bacterium_rate'] * 3.0 * 18.0 / (18 + self.params['m_r_destroys_bacterium_half_sat']))


if __name__ == '__main__':
    unittest.main()
