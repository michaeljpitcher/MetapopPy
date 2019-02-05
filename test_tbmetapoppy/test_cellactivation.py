import unittest
from tbmetapoppy import *


class CellActivationTestCase(unittest.TestCase):

    def setUp(self):
        self.event_mr_t = CellActivation(TBPulmonaryNetwork.MACROPHAGE_RESTING, [TBPulmonaryNetwork.T_CELL_ACTIVATED])
        self.event_mr_b = CellActivation(TBPulmonaryNetwork.MACROPHAGE_RESTING,
                                         TBPulmonaryNetwork.EXTRACELLULAR_BACTERIA)
        self.event_t = CellActivation(TBPulmonaryNetwork.T_CELL_NAIVE,
                                      [TBPulmonaryNetwork.DENDRITIC_CELL_MATURE, TBPulmonaryNetwork.MACROPHAGE_INFECTED])

        self.params = {'m_r_activation_by_t_a_rate': 0.1, 'm_r_activation_by_b_er_b_ed_rate': 0.2,
                       't_n_activation_by_d_m_m_i_rate': 0.3,
                       'm_r_activation_by_t_a_half_sat': 100, 'm_r_activation_by_b_er_b_ed_half_sat': 200,
                       't_n_activation_by_d_m_m_i_half_sat': 300}

        for e in [self.event_mr_t, self.event_mr_b, self.event_t]:
            e.set_parameters(self.params)

        self.network = TBPulmonaryNetwork({TBPulmonaryNetwork.TOPOLOGY: TBPulmonaryNetwork.SINGLE_PATCH})
        self.network.reset()

    def test_rate(self):
        self.assertFalse(self.event_mr_t.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH))
        self.assertFalse(self.event_mr_b.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH))
        self.assertFalse(self.event_t.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH))

        self.network.update_patch(TBPulmonaryNetwork.ALVEOLAR_PATCH,
                                  compartment_changes={TBPulmonaryNetwork.MACROPHAGE_RESTING: 2,
                                                       TBPulmonaryNetwork.T_CELL_NAIVE: 3})
        self.assertFalse(self.event_mr_t.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH))
        self.assertFalse(self.event_mr_b.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH))
        self.assertFalse(self.event_t.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH))

        self.network.update_patch(TBPulmonaryNetwork.ALVEOLAR_PATCH,
                                  compartment_changes={TBPulmonaryNetwork.T_CELL_ACTIVATED: 5})
        self.assertEqual(self.event_mr_t.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH),
                         self.params['m_r_activation_by_t_a_rate'] * 2 *
                         (5.0 / (5 + self.params['m_r_activation_by_t_a_half_sat'])))
        self.assertFalse(self.event_mr_b.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH))
        self.assertFalse(self.event_t.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH))

        self.network.update_patch(TBPulmonaryNetwork.ALVEOLAR_PATCH,
                                  compartment_changes={TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING: 7})
        self.assertEqual(self.event_mr_t.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH),
                         self.params['m_r_activation_by_t_a_rate'] * 2 *
                         (5.0 / (5 + self.params['m_r_activation_by_t_a_half_sat'])))
        self.assertEqual(self.event_mr_b.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH),
                         self.params['m_r_activation_by_b_er_b_ed_rate'] * 2 *
                         (7.0 / (7 + self.params['m_r_activation_by_b_er_b_ed_half_sat'])))
        self.assertFalse(self.event_t.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH))

        self.network.update_patch(TBPulmonaryNetwork.ALVEOLAR_PATCH,
                                  compartment_changes={TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT: 11})
        self.assertEqual(self.event_mr_t.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH),
                         self.params['m_r_activation_by_t_a_rate'] * 2 *
                         (5.0 / (5 + self.params['m_r_activation_by_t_a_half_sat'])))
        self.assertEqual(self.event_mr_b.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH),
                         self.params['m_r_activation_by_b_er_b_ed_rate'] * 2 *
                         (18.0 / (18 + self.params['m_r_activation_by_b_er_b_ed_half_sat'])))
        self.assertFalse(self.event_t.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH))

        self.network.update_patch(TBPulmonaryNetwork.ALVEOLAR_PATCH,
                                  compartment_changes={TBPulmonaryNetwork.DENDRITIC_CELL_MATURE: 13})
        self.assertEqual(self.event_mr_t.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH),
                         self.params['m_r_activation_by_t_a_rate'] * 2 *
                         (5.0 / (5 + self.params['m_r_activation_by_t_a_half_sat'])))
        self.assertEqual(self.event_mr_b.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH),
                         self.params['m_r_activation_by_b_er_b_ed_rate'] * 2 *
                         (18.0 / (18 + self.params['m_r_activation_by_b_er_b_ed_half_sat'])))
        self.assertEqual(self.event_t.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH),
                         self.params['t_n_activation_by_d_m_m_i_rate'] * 3 *
                         (13.0 / (13 + self.params['t_n_activation_by_d_m_m_i_half_sat'])))

        self.network.update_patch(TBPulmonaryNetwork.ALVEOLAR_PATCH,
                                  compartment_changes={TBPulmonaryNetwork.MACROPHAGE_INFECTED: 17})
        self.assertEqual(self.event_mr_t.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH),
                         self.params['m_r_activation_by_t_a_rate'] * 2 *
                         (5.0 / (5 + self.params['m_r_activation_by_t_a_half_sat'])))
        self.assertEqual(self.event_mr_b.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH),
                         self.params['m_r_activation_by_b_er_b_ed_rate'] * 2 *
                         (18.0 / (18 + self.params['m_r_activation_by_b_er_b_ed_half_sat'])))
        self.assertAlmostEqual(self.event_t.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH),
                               self.params['t_n_activation_by_d_m_m_i_rate'] *
                               3 * (30.0 / (30 + self.params['t_n_activation_by_d_m_m_i_half_sat'])))

    def test_perform(self):
        self.network.update_patch(TBPulmonaryNetwork.ALVEOLAR_PATCH,
                                  {TBPulmonaryNetwork.MACROPHAGE_RESTING: 2, TBPulmonaryNetwork.T_CELL_NAIVE: 1})

        self.event_mr_t.perform(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH)
        self.assertEqual(self.network.get_compartment_value(TBPulmonaryNetwork.ALVEOLAR_PATCH,
                                                            TBPulmonaryNetwork.MACROPHAGE_RESTING), 1)
        self.assertEqual(self.network.get_compartment_value(TBPulmonaryNetwork.ALVEOLAR_PATCH,
                                                            TBPulmonaryNetwork.MACROPHAGE_ACTIVATED), 1)
        self.assertEqual(self.network.get_compartment_value(TBPulmonaryNetwork.ALVEOLAR_PATCH,
                                                            TBPulmonaryNetwork.T_CELL_NAIVE), 1)

        self.event_mr_b.perform(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH)
        self.assertEqual(self.network.get_compartment_value(TBPulmonaryNetwork.ALVEOLAR_PATCH,
                                                            TBPulmonaryNetwork.MACROPHAGE_RESTING), 0)
        self.assertEqual(self.network.get_compartment_value(TBPulmonaryNetwork.ALVEOLAR_PATCH,
                                                            TBPulmonaryNetwork.MACROPHAGE_ACTIVATED), 2)
        self.assertEqual(self.network.get_compartment_value(TBPulmonaryNetwork.ALVEOLAR_PATCH,
                                                            TBPulmonaryNetwork.T_CELL_NAIVE), 1)

        self.event_t.perform(self.network, TBPulmonaryNetwork.ALVEOLAR_PATCH)
        self.assertEqual(self.network.get_compartment_value(TBPulmonaryNetwork.ALVEOLAR_PATCH,
                                                            TBPulmonaryNetwork.T_CELL_NAIVE), 0)
        self.assertEqual(self.network.get_compartment_value(TBPulmonaryNetwork.ALVEOLAR_PATCH,
                                                            TBPulmonaryNetwork.T_CELL_ACTIVATED), 1)


if __name__ == '__main__':
    unittest.main()
