import unittest
from tbmetapoppy import *

TEST_CHANGE_RATE_R_TO_D = 'test_change_rate_r_to_d'
TEST_SIGMOID_R_TO_D = 'test_change_sigmoid_r_to_d'
TEST_HALFSAT_R_TO_D = 'test_change_halfsat_r_to_d'
TEST_CHANGE_RATE_D_TO_R = 'test_change_rate_d_to_r'
TEST_SIGMOID_D_TO_R = 'test_change_sigmoid_d_to_r'
TEST_HALFSAT_D_TO_R = 'test_change_halfsat_d_to_r'


class BacteriumChangeStateThroughOxygenTestCase(unittest.TestCase):

    def setUp(self):
        self.event_r_to_d = BacteriumChangeStateThroughOxygen(True, TEST_CHANGE_RATE_R_TO_D, TEST_SIGMOID_R_TO_D,
                                                              TEST_HALFSAT_R_TO_D)
        self.event_d_to_r = BacteriumChangeStateThroughOxygen(False, TEST_CHANGE_RATE_D_TO_R, TEST_SIGMOID_D_TO_R,
                                                              TEST_HALFSAT_D_TO_R)
        self.params = {TEST_CHANGE_RATE_D_TO_R: 0.1, TEST_SIGMOID_D_TO_R: 2, TEST_HALFSAT_D_TO_R: 0.5,
                       TEST_CHANGE_RATE_R_TO_D: 0.2, TEST_SIGMOID_R_TO_D: -3, TEST_HALFSAT_R_TO_D: 0.7}

        self.event_r_to_d.set_parameters(self.params)
        self.event_d_to_r.set_parameters(self.params)

        self.network = TBPulmonaryNetwork({TBPulmonaryNetwork.TOPOLOGY: None})
        self.network.add_node(1)
        self.network.set_patch_type(1, TBPulmonaryNetwork.ALVEOLAR_PATCH)
        self.network.reset()

    def test_rate(self):
        self.network.update_patch(1, attribute_changes={TBPulmonaryNetwork.OXYGEN_TENSION: 1.4})
        self.assertFalse(self.event_d_to_r.calculate_rate_at_patch(self.network, 1))
        self.assertFalse(self.event_r_to_d.calculate_rate_at_patch(self.network, 1))

        self.network.update_patch(1, {TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING: 2, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT: 2})
        r_to_d_1_4 = self.event_r_to_d.calculate_rate_at_patch(self.network, 1)
        self.assertEqual(r_to_d_1_4, self.params[TEST_CHANGE_RATE_R_TO_D] * 2 * (
                    (1.4 ** self.params[TEST_SIGMOID_R_TO_D]) / (
                        1.4 ** self.params[TEST_SIGMOID_R_TO_D] + self.params[TEST_HALFSAT_R_TO_D] ** self.params[
                    TEST_SIGMOID_R_TO_D])))
        d_to_r_1_4 = self.event_d_to_r.calculate_rate_at_patch(self.network, 1)
        self.assertEqual(d_to_r_1_4, self.params[TEST_CHANGE_RATE_D_TO_R] * 2 * (
                (1.4 ** self.params[TEST_SIGMOID_D_TO_R]) / (
                1.4 ** self.params[TEST_SIGMOID_D_TO_R] + self.params[TEST_HALFSAT_D_TO_R] ** self.params[
            TEST_SIGMOID_D_TO_R])))

        self.network.update_patch(1, attribute_changes={TBPulmonaryNetwork.OXYGEN_TENSION: 0.1})
        r_to_d_1_5 = self.event_r_to_d.calculate_rate_at_patch(self.network, 1)
        self.assertEqual(r_to_d_1_5, self.params[TEST_CHANGE_RATE_R_TO_D] * 2 * (
                (1.5 ** self.params[TEST_SIGMOID_R_TO_D]) / (
                1.5 ** self.params[TEST_SIGMOID_R_TO_D] + self.params[TEST_HALFSAT_R_TO_D] ** self.params[
            TEST_SIGMOID_R_TO_D])))
        d_to_r_1_5 = self.event_d_to_r.calculate_rate_at_patch(self.network, 1)
        self.assertEqual(d_to_r_1_5, self.params[TEST_CHANGE_RATE_D_TO_R] * 2 * (
                (1.5 ** self.params[TEST_SIGMOID_D_TO_R]) / (
                1.5 ** self.params[TEST_SIGMOID_D_TO_R] + self.params[TEST_HALFSAT_D_TO_R] ** self.params[
            TEST_SIGMOID_D_TO_R])))

        self.assertTrue(r_to_d_1_4 > r_to_d_1_5) # Oxygen up, so less chance to switch to dormant
        self.assertTrue(d_to_r_1_4 < d_to_r_1_5) # Oxygen down, so more chance to switch to dormant

    def test_perform(self):
        self.network.update_patch(1, {TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING: 10, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT: 10})
        self.event_r_to_d.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING), 9)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT), 11)



if __name__ == '__main__':
    unittest.main()
