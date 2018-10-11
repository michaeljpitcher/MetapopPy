import unittest
from tbmetapoppy import *

TEST_RECRUITMENT_RATE = 'test_recruitment_rate'
TEST_HALF_SAT = 'test_half_sat'

class CellRecruitmentTestCase(unittest.TestCase):

    def setUp(self):
        self.event = CellRecruitment(PulmonaryNetwork.ALVEOLAR_PATCH, TEST_RECRUITMENT_RATE, MACROPHAGE_RESTING)
        self.params = {TEST_RECRUITMENT_RATE: 0.1}
        self.event.set_parameters(self.params)
        self.network = PulmonaryNetwork({PulmonaryNetwork.TOPOLOGY: None}, TB_COMPARTMENTS)
        self.network.add_node(1)
        self.network.set_patch_type(1, PulmonaryNetwork.ALVEOLAR_PATCH)
        self.network.prepare()

    def test_perform(self):
        self.assertEqual(self.network.get_compartment_value(1, MACROPHAGE_RESTING), 0)
        self.event.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, MACROPHAGE_RESTING), 1)


class StandardCellRecruitmentLungTestCase(unittest.TestCase):

    def setUp(self):
        self.event = StandardCellRecruitmentLung(TEST_RECRUITMENT_RATE, MACROPHAGE_RESTING)
        self.params = {TEST_RECRUITMENT_RATE: 0.1}
        self.event.set_parameters(self.params)
        self.network = PulmonaryNetwork({PulmonaryNetwork.TOPOLOGY: None}, TB_COMPARTMENTS)
        self.network.add_node(1)
        self.network.set_patch_type(1, PulmonaryNetwork.ALVEOLAR_PATCH)
        self.network.prepare()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, attribute_changes={PulmonaryNetwork.PERFUSION: 0.7})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.params[TEST_RECRUITMENT_RATE] * 0.7)
        self.network.update_patch(1, attribute_changes={PulmonaryNetwork.PERFUSION: 0.2})
        self.assertAlmostEqual(self.event.calculate_rate_at_patch(self.network, 1),
                               self.params[TEST_RECRUITMENT_RATE] * 0.9)


class EnhancedCellRecruitmentLungTestCase(unittest.TestCase):

    def setUp(self):
        self.event = EnhancedCellRecruitmentLung(TEST_RECRUITMENT_RATE, MACROPHAGE_RESTING,
                                                 BACTERIUM_EXTRACELLULAR_REPLICATING, TEST_HALF_SAT)
        self.params = {TEST_RECRUITMENT_RATE: 0.1, TEST_HALF_SAT: 10}
        self.event.set_parameters(self.params)
        self.network = PulmonaryNetwork({PulmonaryNetwork.TOPOLOGY: None}, TB_COMPARTMENTS)
        self.network.add_node(1)
        self.network.set_patch_type(1, PulmonaryNetwork.ALVEOLAR_PATCH)
        self.network.prepare()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, attribute_changes={PulmonaryNetwork.PERFUSION: 0.7})
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {BACTERIUM_EXTRACELLULAR_REPLICATING: 2})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.params[TEST_RECRUITMENT_RATE] *
                         0.7 * (2.0 / (2 + self.params[TEST_HALF_SAT])))


class StandardCellRecruitmentLymphTestCase(unittest.TestCase):

    def setUp(self):
        self.event = StandardCellRecruitmentLymph(TEST_RECRUITMENT_RATE, T_CELL_NAIVE)
        self.params = {TEST_RECRUITMENT_RATE: 0.1}
        self.event.set_parameters(self.params)
        self.network = PulmonaryNetwork({PulmonaryNetwork.TOPOLOGY: None}, TB_COMPARTMENTS)
        self.network.add_node(1)
        self.network.set_patch_type(1, PulmonaryNetwork.LYMPH_PATCH)
        self.network.prepare()

    def test_rate(self):
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.params[TEST_RECRUITMENT_RATE])


class EnhancedCellRecruitmentLymphTestCase(unittest.TestCase):

    def setUp(self):
        self.event = EnhancedCellRecruitmentLymph(TEST_RECRUITMENT_RATE, MACROPHAGE_RESTING,
                                                  BACTERIUM_EXTRACELLULAR_REPLICATING, TEST_HALF_SAT)
        self.params = {TEST_RECRUITMENT_RATE: 0.1, TEST_HALF_SAT: 10}
        self.event.set_parameters(self.params)
        self.network = PulmonaryNetwork({PulmonaryNetwork.TOPOLOGY: None}, TB_COMPARTMENTS)
        self.network.add_node(1)
        self.network.set_patch_type(1, PulmonaryNetwork.LYMPH_PATCH)
        self.network.prepare()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {BACTERIUM_EXTRACELLULAR_REPLICATING: 2})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.params[TEST_RECRUITMENT_RATE] *
                         (2.0 / (2 + self.params[TEST_HALF_SAT])))


if __name__ == '__main__':
    unittest.main()
