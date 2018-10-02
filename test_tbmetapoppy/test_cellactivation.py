import unittest
from tbmetapoppy import *


class MacrophageActivationTestCase(unittest.TestCase):

    def setUp(self):
        self.event_mr_t = CellActivation(MACROPHAGE_RESTING, [T_CELL_ACTIVATED])
        self.event_mr_b = CellActivation(MACROPHAGE_RESTING, EXTRACELLULAR_BACTERIA)
        self.event_t = CellActivation(T_CELL_NAIVE, [DENDRITIC_CELL_MATURE, MACROPHAGE_INFECTED])

        self.rp = 0.1
        self.hs1 = 100
        self.hs2 = 200
        self.hs3 = 300

        self.event_mr_t.set_parameters(self.rp, self.hs1)
        self.event_mr_b.set_parameters(self.rp, self.hs2)
        self.event_t.set_parameters(self.rp, self.hs3)

        self.network = PulmonaryNetwork({PulmonaryNetwork.TOPOLOGY: None}, TB_COMPARTMENTS)
        self.network.add_node(1)
        self.network.set_patch_type(1, PulmonaryNetwork.ALVEOLAR_PATCH)
        self.network.prepare()

    def test_rate(self):
        self.assertFalse(self.event_mr_t.calculate_rate_at_patch(self.network, 1))
        self.assertFalse(self.event_mr_b.calculate_rate_at_patch(self.network, 1))
        self.assertFalse(self.event_t.calculate_rate_at_patch(self.network, 1))

        self.network.update_patch(1, compartment_changes={MACROPHAGE_RESTING: 2, T_CELL_NAIVE: 3})
        self.assertFalse(self.event_mr_t.calculate_rate_at_patch(self.network, 1))
        self.assertFalse(self.event_mr_b.calculate_rate_at_patch(self.network, 1))
        self.assertFalse(self.event_t.calculate_rate_at_patch(self.network, 1))

        self.network.update_patch(1, compartment_changes={T_CELL_ACTIVATED: 5})
        self.assertEqual(self.event_mr_t.calculate_rate_at_patch(self.network, 1), self.rp * 2 * (5.0 / (5 + self.hs1)))
        self.assertFalse(self.event_mr_b.calculate_rate_at_patch(self.network, 1))
        self.assertFalse(self.event_t.calculate_rate_at_patch(self.network, 1))

        self.network.update_patch(1, compartment_changes={BACTERIUM_EXTRACELLULAR_REPLICATING: 7})
        self.assertEqual(self.event_mr_t.calculate_rate_at_patch(self.network, 1), self.rp * 2 * (5.0 / (5 + self.hs1)))
        self.assertEqual(self.event_mr_b.calculate_rate_at_patch(self.network, 1), self.rp * 2 * (7.0 / (7 + self.hs2)))
        self.assertFalse(self.event_t.calculate_rate_at_patch(self.network, 1))

        self.network.update_patch(1, compartment_changes={BACTERIUM_EXTRACELLULAR_DORMANT: 11})
        self.assertEqual(self.event_mr_t.calculate_rate_at_patch(self.network, 1), self.rp * 2 * (5.0 / (5 + self.hs1)))
        self.assertEqual(self.event_mr_b.calculate_rate_at_patch(self.network, 1), self.rp * 2 * (18.0 / (18 + self.hs2)))
        self.assertFalse(self.event_t.calculate_rate_at_patch(self.network, 1))

        self.network.update_patch(1, compartment_changes={DENDRITIC_CELL_MATURE: 13})
        self.assertEqual(self.event_mr_t.calculate_rate_at_patch(self.network, 1), self.rp * 2 * (5.0 / (5 + self.hs1)))
        self.assertEqual(self.event_mr_b.calculate_rate_at_patch(self.network, 1), self.rp * 2 * (18.0 / (18 + self.hs2)))
        self.assertEqual(self.event_t.calculate_rate_at_patch(self.network, 1), self.rp * 3 * (13.0 / (13 + self.hs3)))

        self.network.update_patch(1, compartment_changes={MACROPHAGE_INFECTED: 17})
        self.assertEqual(self.event_mr_t.calculate_rate_at_patch(self.network, 1), self.rp * 2 * (5.0 / (5 + self.hs1)))
        self.assertEqual(self.event_mr_b.calculate_rate_at_patch(self.network, 1), self.rp * 2 * (18.0 / (18 + self.hs2)))
        self.assertAlmostEqual(self.event_t.calculate_rate_at_patch(self.network, 1), self.rp * 3 * (30.0 / (30 + self.hs3)))

    def test_perform(self):
        self.network.update_patch(1, {MACROPHAGE_RESTING: 2, T_CELL_NAIVE: 1})

        self.event_mr_t.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, MACROPHAGE_RESTING), 1)
        self.assertEqual(self.network.get_compartment_value(1, MACROPHAGE_ACTIVATED), 1)
        self.assertEqual(self.network.get_compartment_value(1, T_CELL_NAIVE), 1)

        self.event_mr_b.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, MACROPHAGE_RESTING), 0)
        self.assertEqual(self.network.get_compartment_value(1, MACROPHAGE_ACTIVATED), 2)
        self.assertEqual(self.network.get_compartment_value(1, T_CELL_NAIVE), 1)

        self.event_t.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, T_CELL_NAIVE), 0)
        self.assertEqual(self.network.get_compartment_value(1, T_CELL_ACTIVATED), 1)


if __name__ == '__main__':
    unittest.main()
