import unittest
from tbmetapoppy import *


class BacteriumChangeStateThroughOxygenTestCase(unittest.TestCase):

    def setUp(self):
        self.event_r_to_d = BacteriumChangeStateThroughOxygen(True)
        self.event_d_to_r = BacteriumChangeStateThroughOxygen(False)

        self.params = {BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_REPLICATING + RATE: 0.1,
                       BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_REPLICATING + SIGMOID: 0.2,
                       BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_REPLICATING + HALF_SAT: 0.3,
                       BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_DORMANT + RATE: 7.4,
                       BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_DORMANT + SIGMOID: -7.5,
                       BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_DORMANT + HALF_SAT: 7.6}

        self.event_r_to_d.set_parameters(self.params)
        self.event_d_to_r.set_parameters(self.params)

        self.network = TBPulmonaryEnvironment({TBPulmonaryEnvironment.TOPOLOGY: None})
        self.network.add_node(1)
        self.network.set_patch_type(1, TBPulmonaryEnvironment.ALVEOLAR_PATCH)
        self.network.reset()

    def test_rate(self):
        self.network.update_patch(1, attribute_changes={TBPulmonaryEnvironment.OXYGEN_TENSION: 1.4})
        self.assertFalse(self.event_d_to_r.calculate_rate_at_patch(self.network, 1))
        self.assertFalse(self.event_r_to_d.calculate_rate_at_patch(self.network, 1))

        self.network.update_patch(1, {TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_REPLICATING: 2, TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT: 2})
        r_to_d_1_4 = self.event_r_to_d.calculate_rate_at_patch(self.network, 1)
        self.assertEqual(r_to_d_1_4, self.params[BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_DORMANT + RATE] * 2 * (
                    (1.4 ** self.params[BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_DORMANT + SIGMOID]) / (
                        1.4 ** self.params[BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_DORMANT + SIGMOID] + self.params[BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_DORMANT + HALF_SAT] ** self.params[
            BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_DORMANT + SIGMOID])))
        d_to_r_1_4 = self.event_d_to_r.calculate_rate_at_patch(self.network, 1)
        self.assertEqual(d_to_r_1_4, self.params[BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_REPLICATING + RATE] * 2 * (
                (1.4 ** self.params[BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_REPLICATING + SIGMOID]) / (
                1.4 ** self.params[BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_REPLICATING + SIGMOID] + self.params[BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_REPLICATING + HALF_SAT] ** self.params[
            BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_REPLICATING + SIGMOID])))

        self.network.update_patch(1, attribute_changes={TBPulmonaryEnvironment.OXYGEN_TENSION: 0.1})
        r_to_d_1_5 = self.event_r_to_d.calculate_rate_at_patch(self.network, 1)
        self.assertEqual(r_to_d_1_5, self.params[BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_DORMANT + RATE] * 2 * (
                (1.5 ** self.params[BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_DORMANT + SIGMOID]) / (
                1.5 ** self.params[BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_DORMANT + SIGMOID] + self.params[BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_DORMANT + HALF_SAT] ** self.params[
            BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_DORMANT + SIGMOID])))
        d_to_r_1_5 = self.event_d_to_r.calculate_rate_at_patch(self.network, 1)
        self.assertEqual(d_to_r_1_5, self.params[BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_REPLICATING + RATE] * 2 * (
                (1.5 ** self.params[BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_REPLICATING + SIGMOID]) / (
                1.5 ** self.params[BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_REPLICATING + SIGMOID] + self.params[BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_REPLICATING + HALF_SAT] ** self.params[
            BacteriumChangeStateThroughOxygen.BACTERIUM_CHANGE_TO_REPLICATING + SIGMOID])))

        self.assertTrue(r_to_d_1_4 > r_to_d_1_5) # Oxygen up, so less chance to switch to dormant
        self.assertTrue(d_to_r_1_4 < d_to_r_1_5) # Oxygen down, so more chance to switch to dormant

    def test_perform(self):
        self.network.update_patch(1, {TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_REPLICATING: 10, TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT: 10})
        self.event_r_to_d.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_REPLICATING), 9)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT), 11)



if __name__ == '__main__':
    unittest.main()
