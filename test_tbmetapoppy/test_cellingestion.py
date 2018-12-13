import unittest
from tbmetapoppy import *

TEST_RATE_MR_INGESTION = 'test_rate_mr_ingestion'
TEST_HALFSAT_MR_INGESTION = 'test_halfsat_mr_ingestion'
TEST_PROB_MR_INFECTION = 'test_prob_mr_infection'

TEST_RATE_MA_INGESTION = 'test_rate_ma_ingestion'
TEST_HALFSAT_MA_INGESTION = 'test_halfsat_ma_ingestion'
TEST_PROB_MA_INFECTION = 'test_prob_ma_infection'

TEST_RATE_DC_INGESTION = 'test_rate_mac_ingestion'
TEST_HALFSAT_DC_INGESTION = 'test_halfsat_mac_ingestion'
TEST_PROB_DC_INFECTION = 'test_prob_mac_infection'


class CellIngestBacteriumTestCase(unittest.TestCase):

    def setUp(self):
        self.event_mr = CellIngestBacterium(TEST_RATE_MR_INGESTION, TBPulmonaryNetwork.MACROPHAGE_RESTING, TEST_HALFSAT_MR_INGESTION,
                                             TEST_PROB_MR_INFECTION)
        self.event_ma = CellIngestBacterium(TEST_RATE_MA_INGESTION, TBPulmonaryNetwork.MACROPHAGE_ACTIVATED, TEST_HALFSAT_MA_INGESTION,
                                             TEST_PROB_MA_INFECTION)
        self.event_dc = CellIngestBacterium(TEST_RATE_DC_INGESTION, TBPulmonaryNetwork.DENDRITIC_CELL_IMMATURE, TEST_HALFSAT_DC_INGESTION,
                                             TEST_PROB_DC_INFECTION)

        self.params = {TEST_RATE_MR_INGESTION: 0.1, TEST_HALFSAT_MR_INGESTION: 101, TEST_PROB_MR_INFECTION: 0.7,
                       TEST_RATE_MA_INGESTION: 0.2, TEST_HALFSAT_MA_INGESTION: 99, TEST_PROB_MA_INFECTION: 0.0,
                       TEST_RATE_DC_INGESTION: 0.1, TEST_HALFSAT_DC_INGESTION: 101, TEST_PROB_DC_INFECTION: 1.0}

        self.events = [self.event_mr, self.event_ma, self.event_dc]
        for e in self.events:
            e.set_parameters(self.params)
        self.network = TBPulmonaryNetwork({TBPulmonaryNetwork.TOPOLOGY: None})
        self.network.add_node(1)
        self.network.set_patch_type(1, TBPulmonaryNetwork.ALVEOLAR_PATCH)
        self.network.reset()

    def test_rate(self):
        for e in self.events:
            self.assertFalse(e.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {TBPulmonaryNetwork.MACROPHAGE_RESTING: 3, TBPulmonaryNetwork.MACROPHAGE_ACTIVATED: 5, TBPulmonaryNetwork.DENDRITIC_CELL_IMMATURE: 7})
        for e in self.events:
            self.assertFalse(e.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING: 11})
        self.assertEqual(self.event_mr.calculate_rate_at_patch(self.network, 1),
                         self.params[TEST_RATE_MR_INGESTION] * 3 * (
                                 11.0 / (11.0 + self.params[TEST_HALFSAT_MR_INGESTION])))
        self.assertEqual(self.event_ma.calculate_rate_at_patch(self.network, 1),
                         self.params[TEST_RATE_MA_INGESTION] * 5 * (
                                     11.0 / (11.0 + self.params[TEST_HALFSAT_MA_INGESTION])))
        self.assertEqual(self.event_dc.calculate_rate_at_patch(self.network, 1),
                         self.params[TEST_RATE_DC_INGESTION] * 7 * (
                                     11.0 / (11.0 + self.params[TEST_HALFSAT_DC_INGESTION])))
        self.network.update_patch(1, {TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT: 13})
        self.assertEqual(self.event_mr.calculate_rate_at_patch(self.network, 1),
                         self.params[TEST_RATE_MR_INGESTION] * 3 * (
                                 24.0 / (24.0 + self.params[TEST_HALFSAT_MR_INGESTION])))
        self.assertEqual(self.event_ma.calculate_rate_at_patch(self.network, 1),
                         self.params[TEST_RATE_MA_INGESTION] * 5 * (
                                 24.0 / (24.0 + self.params[TEST_HALFSAT_MA_INGESTION])))
        self.assertEqual(self.event_dc.calculate_rate_at_patch(self.network, 1),
                         self.params[TEST_RATE_DC_INGESTION] * 7 * (
                                 24.0 / (24.0 + self.params[TEST_HALFSAT_DC_INGESTION])))

    def test_perform(self):
        self.network.update_patch(1, {TBPulmonaryNetwork.MACROPHAGE_RESTING: 10, TBPulmonaryNetwork.MACROPHAGE_ACTIVATED:10, TBPulmonaryNetwork.DENDRITIC_CELL_IMMATURE:10,
                                      TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING:15, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT: 15})
        # Never infected
        for n in range(10):
            self.event_ma.perform(self.network, 1)
            self.assertEqual(self.network.get_compartment_value(1, [TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT,
                                                                    TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING]), 30 - n - 1)
        # Always infected
        for n in range(10):
            self.event_dc.perform(self.network, 1)
            self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.DENDRITIC_CELL_IMMATURE), 10-n-1)
            self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_DENDRITIC), n+1)

        # Either
        for n in range(10):
            self.event_mr.perform(self.network, 1)
            self.assertEqual(self.network.get_compartment_value(1, [TBPulmonaryNetwork.MACROPHAGE_RESTING, TBPulmonaryNetwork.MACROPHAGE_INFECTED]), 10)
            self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.MACROPHAGE_INFECTED),
                             self.network.get_compartment_value(1, TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE))



if __name__ == '__main__':
    unittest.main()
