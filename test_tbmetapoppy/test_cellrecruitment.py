import unittest
from tbmetapoppy import *

class StandardCellRecruitmentLungTestCase(unittest.TestCase):

    def setUp(self):
        self.event = StandardCellRecruitmentLung(TBPulmonaryEnvironment.MACROPHAGE_RESTING)
        self.params = {'m_r_standard_recruitment_alveolar_patch_rate': 0.1}
        self.event.set_parameters(self.params)
        self.network = TBPulmonaryEnvironment({TBPulmonaryEnvironment.TOPOLOGY: TBPulmonaryEnvironment.SINGLE_PATCH})
        self.network.reset()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, TBPulmonaryEnvironment.ALVEOLAR_PATCH))
        self.network.update_patch(TBPulmonaryEnvironment.ALVEOLAR_PATCH, attribute_changes={TBPulmonaryEnvironment.PERFUSION: 0.7})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, TBPulmonaryEnvironment.ALVEOLAR_PATCH),
                         self.params['m_r_standard_recruitment_alveolar_patch_rate'] * 0.7)
        self.network.update_patch(TBPulmonaryEnvironment.ALVEOLAR_PATCH, attribute_changes={TBPulmonaryEnvironment.PERFUSION: 0.2})
        self.assertAlmostEqual(self.event.calculate_rate_at_patch(self.network, TBPulmonaryEnvironment.ALVEOLAR_PATCH),
                               self.params['m_r_standard_recruitment_alveolar_patch_rate'] * 0.9)

    def test_perform(self):
        self.assertEqual(self.network.get_compartment_value(TBPulmonaryEnvironment.ALVEOLAR_PATCH, TBPulmonaryEnvironment.MACROPHAGE_RESTING), 0)
        self.event.perform(self.network, TBPulmonaryEnvironment.ALVEOLAR_PATCH)
        self.assertEqual(self.network.get_compartment_value(TBPulmonaryEnvironment.ALVEOLAR_PATCH, TBPulmonaryEnvironment.MACROPHAGE_RESTING), 1)


class EnhancedCellRecruitmentLungTestCase(unittest.TestCase):

    def setUp(self):
        self.event = EnhancedCellRecruitmentLung(TBPulmonaryEnvironment.MACROPHAGE_RESTING)
        self.params = {'m_r_enhanced_recruitment_alveolar_patch_rate': 0.1,
                       'm_r_enhanced_recruitment_alveolar_patch_half_sat': 10,
                       'macrophage_infected_to_activated_chemokine_weight': 0.2}
        self.event.set_parameters(self.params)
        self.network = TBPulmonaryEnvironment({TBPulmonaryEnvironment.TOPOLOGY: TBPulmonaryEnvironment.SINGLE_PATCH})
        self.network.reset()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, TBPulmonaryEnvironment.ALVEOLAR_PATCH))
        self.network.update_patch(TBPulmonaryEnvironment.ALVEOLAR_PATCH, attribute_changes={TBPulmonaryEnvironment.PERFUSION: 0.7})
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, TBPulmonaryEnvironment.ALVEOLAR_PATCH))
        self.network.update_patch(TBPulmonaryEnvironment.ALVEOLAR_PATCH, {TBPulmonaryEnvironment.MACROPHAGE_INFECTED: 2,
                                                                          TBPulmonaryEnvironment.MACROPHAGE_ACTIVATED: 3})
        self.assertAlmostEqual(self.event.calculate_rate_at_patch(self.network, TBPulmonaryEnvironment.ALVEOLAR_PATCH),
                               0.1 * 0.7 * (3.0 + (0.2 * 2)) / (3.0 + (0.2 * 2) + 10))


class StandardCellRecruitmentLymphTestCase(unittest.TestCase):

    def setUp(self):
        self.event = StandardCellRecruitmentLymph(TBPulmonaryEnvironment.T_CELL_NAIVE)
        self.params = {'t_n_standard_recruitment_lymph_patch_rate': 0.1}
        self.event.set_parameters(self.params)
        self.network = TBPulmonaryEnvironment({TBPulmonaryEnvironment.TOPOLOGY: TBPulmonaryEnvironment.SINGLE_PATCH})
        self.network.reset()

    def test_rate(self):
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, TBPulmonaryEnvironment.LYMPH_PATCH),
                         self.params['t_n_standard_recruitment_lymph_patch_rate'])


class EnhancedCellRecruitmentLymphTestCase(unittest.TestCase):

    def setUp(self):
        self.event = EnhancedCellRecruitmentLymph(TBPulmonaryEnvironment.MACROPHAGE_RESTING)
        self.params = {'m_r_enhanced_recruitment_lymph_patch_rate': 0.1,
                       'm_r_enhanced_recruitment_lymph_patch_half_sat': 10,
                       'macrophage_infected_to_activated_chemokine_weight': 0.3}
        self.event.set_parameters(self.params)
        self.network = TBPulmonaryEnvironment({TBPulmonaryEnvironment.TOPOLOGY: TBPulmonaryEnvironment.SINGLE_PATCH})
        self.network.reset()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, TBPulmonaryEnvironment.LYMPH_PATCH))
        self.network.update_patch(TBPulmonaryEnvironment.LYMPH_PATCH, {TBPulmonaryEnvironment.MACROPHAGE_ACTIVATED: 2,
                                                                       TBPulmonaryEnvironment.MACROPHAGE_INFECTED:3})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, TBPulmonaryEnvironment.LYMPH_PATCH),
                         self.params['m_r_enhanced_recruitment_lymph_patch_rate'] *
                         ((2.0 + 0.3*3) / (2 + 0.3*3 + self.params['m_r_enhanced_recruitment_lymph_patch_half_sat'])))


if __name__ == '__main__':
    unittest.main()
