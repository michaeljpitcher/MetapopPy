import unittest
from tbmetapoppy import *

TEST_RECRUITMENT_RATE = 'test_recruitment_rate'
TEST_HALF_SAT = 'test_half_sat'
TEST_WEIGHT = 'test_weight'


class CellRecruitmentTestCase(unittest.TestCase):

    def setUp(self):
        self.event = CellRecruitment(TBPulmonaryNetwork.ALVEOLAR_PATCH, TEST_RECRUITMENT_RATE, TBPulmonaryNetwork.MACROPHAGE_RESTING)
        self.params = {TEST_RECRUITMENT_RATE: 0.1}
        self.event.set_parameters(self.params)
        self.network = TBPulmonaryNetwork({TBPulmonaryNetwork.TOPOLOGY: None})
        self.network.add_node(1)
        self.network.set_patch_type(1, TBPulmonaryNetwork.ALVEOLAR_PATCH)
        self.network.reset()

    def test_perform(self):
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.MACROPHAGE_RESTING), 0)
        self.event.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.MACROPHAGE_RESTING), 1)


class StandardCellRecruitmentLungTestCase(unittest.TestCase):

    def setUp(self):
        self.event = StandardCellRecruitmentLung(TEST_RECRUITMENT_RATE, TBPulmonaryNetwork.MACROPHAGE_RESTING)
        self.params = {TEST_RECRUITMENT_RATE: 0.1}
        self.event.set_parameters(self.params)
        self.network = TBPulmonaryNetwork({TBPulmonaryNetwork.TOPOLOGY: None})
        self.network.add_node(1)
        self.network.set_patch_type(1, TBPulmonaryNetwork.ALVEOLAR_PATCH)
        self.network.reset()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, attribute_changes={TBPulmonaryNetwork.PERFUSION: 0.7})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.params[TEST_RECRUITMENT_RATE] * 0.7)
        self.network.update_patch(1, attribute_changes={TBPulmonaryNetwork.PERFUSION: 0.2})
        self.assertAlmostEqual(self.event.calculate_rate_at_patch(self.network, 1),
                               self.params[TEST_RECRUITMENT_RATE] * 0.9)


class EnhancedCellRecruitmentLungTestCase(unittest.TestCase):

    def setUp(self):
        self.event = EnhancedCellRecruitmentLung(TEST_RECRUITMENT_RATE, TBPulmonaryNetwork.MACROPHAGE_RESTING,
                                                 TEST_HALF_SAT, TEST_WEIGHT)
        self.params = {TEST_RECRUITMENT_RATE: 0.1, TEST_HALF_SAT: 10, TEST_WEIGHT: 0.2}
        self.event.set_parameters(self.params)
        self.network = TBPulmonaryNetwork({TBPulmonaryNetwork.TOPOLOGY: None})
        self.network.add_node(1)
        self.network.set_patch_type(1, TBPulmonaryNetwork.ALVEOLAR_PATCH)
        self.network.reset()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, attribute_changes={TBPulmonaryNetwork.PERFUSION: 0.7})
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {TBPulmonaryNetwork.MACROPHAGE_INFECTED: 2, TBPulmonaryNetwork.MACROPHAGE_ACTIVATED: 3})
        self.assertAlmostEqual(self.event.calculate_rate_at_patch(self.network, 1),
                               0.1 * 0.7 * (3.0 + (0.2 * 2)) / (3.0 + (0.2 * 2) + 10))


class StandardCellRecruitmentLymphTestCase(unittest.TestCase):

    def setUp(self):
        self.event = StandardCellRecruitmentLymph(TEST_RECRUITMENT_RATE, TBPulmonaryNetwork.T_CELL_NAIVE)
        self.params = {TEST_RECRUITMENT_RATE: 0.1}
        self.event.set_parameters(self.params)
        self.network = TBPulmonaryNetwork({TBPulmonaryNetwork.TOPOLOGY: None})
        self.network.add_node(1)
        self.network.set_patch_type(1, TBPulmonaryNetwork.LYMPH_PATCH)
        self.network.reset()

    def test_rate(self):
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.params[TEST_RECRUITMENT_RATE])


class EnhancedCellRecruitmentLymphTestCase(unittest.TestCase):

    def setUp(self):
        self.event = EnhancedCellRecruitmentLymph(TEST_RECRUITMENT_RATE, TBPulmonaryNetwork.MACROPHAGE_RESTING,
                                                  TEST_HALF_SAT, TEST_WEIGHT)
        self.params = {TEST_RECRUITMENT_RATE: 0.1, TEST_HALF_SAT: 10, TEST_WEIGHT: 0.3}
        self.event.set_parameters(self.params)
        self.network = TBPulmonaryNetwork({TBPulmonaryNetwork.TOPOLOGY: None})
        self.network.add_node(1)
        self.network.set_patch_type(1, TBPulmonaryNetwork.LYMPH_PATCH)
        self.network.reset()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {TBPulmonaryNetwork.MACROPHAGE_ACTIVATED: 2, TBPulmonaryNetwork.MACROPHAGE_INFECTED:3})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.params[TEST_RECRUITMENT_RATE] *
                         ((2.0 + 0.3*3) / (2 + 0.3*3 + self.params[TEST_HALF_SAT])))


if __name__ == '__main__':
    unittest.main()
