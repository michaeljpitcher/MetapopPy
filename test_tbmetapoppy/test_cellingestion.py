import unittest
from tbmetapoppy import *



class CellIngestBacteriumTestCase(unittest.TestCase):

    def setUp(self):
        self.event_mr = CellIngestBacterium(TBPulmonaryNetwork.MACROPHAGE_RESTING)
        self.event_ma = CellIngestBacterium(TBPulmonaryNetwork.MACROPHAGE_ACTIVATED)
        self.event_dc = CellIngestBacterium(TBPulmonaryNetwork.DENDRITIC_CELL_IMMATURE)

        self.params = {'m_r_ingest_bacterium_rate': 0.1, 'm_r_ingest_bacterium_half_sat': 101, 'm_r_infection_probability': 0.7,
                       'm_a_ingest_bacterium_rate': 0.2, 'm_a_ingest_bacterium_half_sat': 99, 'm_a_infection_probability': 0.0,
                       'd_i_ingest_bacterium_rate': 0.1, 'd_i_ingest_bacterium_half_sat': 101, 'd_i_infection_probability': 1.0}

        self.events = [self.event_mr, self.event_ma, self.event_dc]
        for e in self.events:
            e.set_parameters(self.params)
        self.network = TBPulmonaryNetwork({TBPulmonaryNetwork.TOPOLOGY: TBPulmonaryNetwork.SINGLE_PATCH})
        self.network.reset()

    def test_rate(self):
        for e in self.events:
            self.assertFalse(e.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH))
        self.network.update_patch(TBPulmonaryNetwork.ALVEOLAR_PATCH, {TBPulmonaryNetwork.MACROPHAGE_RESTING: 3, TBPulmonaryNetwork.MACROPHAGE_ACTIVATED: 5, TBPulmonaryNetwork.DENDRITIC_CELL_IMMATURE: 7})
        for e in self.events:
            self.assertFalse(e.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH))
        self.network.update_patch(TBPulmonaryNetwork.ALVEOLAR_PATCH, {TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING: 11})
        self.assertEqual(self.event_mr.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH),
                         self.params['m_r_ingest_bacterium_rate'] * 3 * (
                                 11.0 / (11.0 + self.params['m_r_ingest_bacterium_half_sat'])))
        self.assertEqual(self.event_ma.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH),
                         self.params['m_a_ingest_bacterium_rate'] * 5 * (
                                     11.0 / (11.0 + self.params['m_a_ingest_bacterium_half_sat'])))
        self.assertEqual(self.event_dc.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH),
                         self.params['d_i_ingest_bacterium_rate'] * 7 * (
                                     11.0 / (11.0 + self.params['d_i_ingest_bacterium_half_sat'])))
        self.network.update_patch(TBPulmonaryNetwork.ALVEOLAR_PATCH, {TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT: 13})
        self.assertEqual(self.event_mr.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH),
                         self.params['m_r_ingest_bacterium_rate'] * 3 * (
                                 24.0 / (24.0 + self.params['m_r_ingest_bacterium_half_sat'])))
        self.assertEqual(self.event_ma.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH),
                         self.params['m_a_ingest_bacterium_rate'] * 5 * (
                                 24.0 / (24.0 + self.params['m_a_ingest_bacterium_half_sat'])))
        self.assertEqual(self.event_dc.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH),
                         self.params['d_i_ingest_bacterium_rate'] * 7 * (
                                 24.0 / (24.0 + self.params['d_i_ingest_bacterium_half_sat'])))

    def test_perform(self):
        self.network.update_patch(TBPulmonaryNetwork.ALVEOLAR_PATCH, {TBPulmonaryNetwork.MACROPHAGE_RESTING: 10, TBPulmonaryNetwork.MACROPHAGE_ACTIVATED:10, TBPulmonaryNetwork.DENDRITIC_CELL_IMMATURE:10,
                                      TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING:15, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT: 15})
        # Never infected
        for n in range(10):
            self.event_ma.perform(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH)
            self.assertEqual(self.network.get_compartment_value(TBPulmonaryNetwork.ALVEOLAR_PATCH, [TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT,
                                                                    TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING]), 30 - n - 1)
        # Always infected
        for n in range(10):
            self.event_dc.perform(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH)
            self.assertEqual(self.network.get_compartment_value(TBPulmonaryNetwork.ALVEOLAR_PATCH, TBPulmonaryNetwork.DENDRITIC_CELL_IMMATURE), 10-n-1)
            self.assertEqual(self.network.get_compartment_value(TBPulmonaryNetwork.ALVEOLAR_PATCH, TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_DENDRITIC), n+1)

        # Either
        for n in range(10):
            self.event_mr.perform(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH)
            self.assertEqual(self.network.get_compartment_value(TBPulmonaryNetwork.ALVEOLAR_PATCH, [TBPulmonaryNetwork.MACROPHAGE_RESTING, TBPulmonaryNetwork.MACROPHAGE_INFECTED]), 10)
            self.assertEqual(self.network.get_compartment_value(TBPulmonaryNetwork.ALVEOLAR_PATCH, TBPulmonaryNetwork.MACROPHAGE_INFECTED),
                             self.network.get_compartment_value(TBPulmonaryNetwork.ALVEOLAR_PATCH, TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE))



if __name__ == '__main__':
    unittest.main()
