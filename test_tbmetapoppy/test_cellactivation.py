import unittest
from tbmetapoppy import *

TEST_RATE_MR_TA = 'test_rate_MR_TA'
TEST_RATE_MR_BE = 'test_rate_MR_BE'
TEST_RATE_TA_APC = 'test_rate_MR_APC'
TEST_HALFSAT_MR_TA = 'test_halfsat_MR_TA'
TEST_HALFSAT_MR_BE = 'test_halfsat_MR_BE'
TEST_HALFSAT_TA_APC = 'test_halfsat_MR_APC'


class CellActivationTestCase(unittest.TestCase):

    def setUp(self):
        self.event_mr_t = CellActivation(TEST_RATE_MR_TA, TEST_HALFSAT_MR_TA, TBPulmonaryNetwork.MACROPHAGE_RESTING, [TBPulmonaryNetwork.T_CELL_ACTIVATED])
        self.event_mr_b = CellActivation(TEST_RATE_MR_BE, TEST_HALFSAT_MR_BE, TBPulmonaryNetwork.MACROPHAGE_RESTING,
                                         TBPulmonaryNetwork.EXTRACELLULAR_BACTERIA)
        self.event_t = CellActivation(TEST_RATE_TA_APC, TEST_HALFSAT_TA_APC, TBPulmonaryNetwork.T_CELL_NAIVE,
                                      [TBPulmonaryNetwork.DENDRITIC_CELL_MATURE, TBPulmonaryNetwork.MACROPHAGE_INFECTED])

        self.params = {TEST_RATE_MR_TA: 0.1, TEST_RATE_MR_BE: 0.2, TEST_RATE_TA_APC: 0.3, TEST_HALFSAT_MR_TA: 100,
                       TEST_HALFSAT_MR_BE: 200, TEST_HALFSAT_TA_APC: 300}

        for e in [self.event_mr_t, self.event_mr_b, self.event_t]:
            e.set_parameters(self.params)

        self.network = TBPulmonaryNetwork({TBPulmonaryNetwork.TOPOLOGY: None})
        self.network.add_node(1)
        self.network.set_patch_type(1, TBPulmonaryNetwork.ALVEOLAR_PATCH)
        self.network.reset()

    def test_rate(self):
        self.assertFalse(self.event_mr_t.calculate_rate_at_patch(self.network, 1))
        self.assertFalse(self.event_mr_b.calculate_rate_at_patch(self.network, 1))
        self.assertFalse(self.event_t.calculate_rate_at_patch(self.network, 1))

        self.network.update_patch(1, compartment_changes={TBPulmonaryNetwork.MACROPHAGE_RESTING: 2, TBPulmonaryNetwork.T_CELL_NAIVE: 3})
        self.assertFalse(self.event_mr_t.calculate_rate_at_patch(self.network, 1))
        self.assertFalse(self.event_mr_b.calculate_rate_at_patch(self.network, 1))
        self.assertFalse(self.event_t.calculate_rate_at_patch(self.network, 1))

        self.network.update_patch(1, compartment_changes={TBPulmonaryNetwork.T_CELL_ACTIVATED: 5})
        self.assertEqual(self.event_mr_t.calculate_rate_at_patch(self.network, 1), self.params[TEST_RATE_MR_TA] * 2 *
                         (5.0 / (5 + self.params[TEST_HALFSAT_MR_TA])))
        self.assertFalse(self.event_mr_b.calculate_rate_at_patch(self.network, 1))
        self.assertFalse(self.event_t.calculate_rate_at_patch(self.network, 1))

        self.network.update_patch(1, compartment_changes={TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING: 7})
        self.assertEqual(self.event_mr_t.calculate_rate_at_patch(self.network, 1), self.params[TEST_RATE_MR_TA] * 2 *
                         (5.0 / (5 + self.params[TEST_HALFSAT_MR_TA])))
        self.assertEqual(self.event_mr_b.calculate_rate_at_patch(self.network, 1), self.params[TEST_RATE_MR_BE] * 2 *
                         (7.0 / (7 + self.params[TEST_HALFSAT_MR_BE])))
        self.assertFalse(self.event_t.calculate_rate_at_patch(self.network, 1))

        self.network.update_patch(1, compartment_changes={TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT: 11})
        self.assertEqual(self.event_mr_t.calculate_rate_at_patch(self.network, 1), self.params[TEST_RATE_MR_TA] * 2 *
                         (5.0 / (5 + self.params[TEST_HALFSAT_MR_TA])))
        self.assertEqual(self.event_mr_b.calculate_rate_at_patch(self.network, 1), self.params[TEST_RATE_MR_BE] * 2 *
                         (18.0 / (18 + self.params[TEST_HALFSAT_MR_BE])))
        self.assertFalse(self.event_t.calculate_rate_at_patch(self.network, 1))

        self.network.update_patch(1, compartment_changes={TBPulmonaryNetwork.DENDRITIC_CELL_MATURE: 13})
        self.assertEqual(self.event_mr_t.calculate_rate_at_patch(self.network, 1), self.params[TEST_RATE_MR_TA] * 2 *
                         (5.0 / (5 + self.params[TEST_HALFSAT_MR_TA])))
        self.assertEqual(self.event_mr_b.calculate_rate_at_patch(self.network, 1), self.params[TEST_RATE_MR_BE] * 2 *
                         (18.0 / (18 + self.params[TEST_HALFSAT_MR_BE])))
        self.assertEqual(self.event_t.calculate_rate_at_patch(self.network, 1), self.params[TEST_RATE_TA_APC] * 3 *
                         (13.0 / (13 + self.params[TEST_HALFSAT_TA_APC])))

        self.network.update_patch(1, compartment_changes={TBPulmonaryNetwork.MACROPHAGE_INFECTED: 17})
        self.assertEqual(self.event_mr_t.calculate_rate_at_patch(self.network, 1), self.params[TEST_RATE_MR_TA] * 2 *
                         (5.0 / (5 + self.params[TEST_HALFSAT_MR_TA])))
        self.assertEqual(self.event_mr_b.calculate_rate_at_patch(self.network, 1), self.params[TEST_RATE_MR_BE] * 2 *
                         (18.0 / (18 + self.params[TEST_HALFSAT_MR_BE])))
        self.assertAlmostEqual(self.event_t.calculate_rate_at_patch(self.network, 1), self.params[TEST_RATE_TA_APC] *
                               3 * (30.0 / (30 + self.params[TEST_HALFSAT_TA_APC])))

    def test_perform(self):
        self.network.update_patch(1, {TBPulmonaryNetwork.MACROPHAGE_RESTING: 2, TBPulmonaryNetwork.T_CELL_NAIVE: 1})

        self.event_mr_t.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.MACROPHAGE_RESTING), 1)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.MACROPHAGE_ACTIVATED), 1)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.T_CELL_NAIVE), 1)

        self.event_mr_b.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.MACROPHAGE_RESTING), 0)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.MACROPHAGE_ACTIVATED), 2)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.T_CELL_NAIVE), 1)

        self.event_t.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.T_CELL_NAIVE), 0)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.T_CELL_ACTIVATED), 1)


if __name__ == '__main__':
    unittest.main()
