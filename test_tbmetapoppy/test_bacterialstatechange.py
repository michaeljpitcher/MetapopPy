import unittest
from tbmetapoppy import *


class BacteriumBecomesDormantTestCase(unittest.TestCase):

    def setUp(self):
        self.event = BacteriumBecomesDormant()
        self.rp = 0.1
        self.sig = -2
        self.hs = 0.5
        self.event.set_parameters(self.rp, self.sig, self.hs)
        self.network = PulmonaryNetwork({PulmonaryNetwork.TOPOLOGY: None}, TB_COMPARTMENTS)
        self.network.add_node(1)
        self.network.set_patch_type(1, PulmonaryNetwork.ALVEOLAR_PATCH)
        self.network.prepare()

    def test_rate(self):
        self.network.update_patch(1,attribute_changes={PulmonaryNetwork.OXYGEN_TENSION: 1.4})
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {BACTERIUM_EXTRACELLULAR_REPLICATING: 2})
        o2_at_1_4 = self.event.calculate_rate_at_patch(self.network, 1)
        self.assertEqual(o2_at_1_4, self.rp * 2 * ((1.4**self.sig)/(1.4**self.sig + self.hs**self.sig)))

        self.network.update_patch(1, attribute_changes={PulmonaryNetwork.OXYGEN_TENSION: -0.4})
        o2_at_1_0 = self.event.calculate_rate_at_patch(self.network, 1)
        self.assertEqual(o2_at_1_0, self.rp * 2 * ((1.0 ** self.sig) / (1.0 ** self.sig + self.hs ** self.sig)))

        # Drop in oxygen increases the rate of change
        self.assertTrue(o2_at_1_0 > o2_at_1_4)

    def test_perform(self):
        self.network.update_patch(1, {BACTERIUM_EXTRACELLULAR_REPLICATING: 1})
        self.event.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_EXTRACELLULAR_REPLICATING), 0)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_EXTRACELLULAR_DORMANT), 1)


class BacteriumBecomesReplicatingTestCase(unittest.TestCase):

    def setUp(self):
        self.event = BacteriumBecomesReplicating()
        self.rp = 0.1
        self.sig = 2
        self.hs = 0.5
        self.event.set_parameters(self.rp, self.sig, self.hs)
        self.network = PulmonaryNetwork({PulmonaryNetwork.TOPOLOGY: None}, TB_COMPARTMENTS)
        self.network.add_node(1)
        self.network.set_patch_type(1, PulmonaryNetwork.ALVEOLAR_PATCH)
        self.network.prepare()

    def test_rate(self):
        self.network.update_patch(1,attribute_changes={PulmonaryNetwork.OXYGEN_TENSION: 0.4})
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {BACTERIUM_EXTRACELLULAR_DORMANT: 2})
        o2_at_0_4 = self.event.calculate_rate_at_patch(self.network, 1)
        self.assertEqual(o2_at_0_4, self.rp * 2 * ((0.4**self.sig)/(0.4**self.sig + self.hs**self.sig)))

        self.network.update_patch(1, attribute_changes={PulmonaryNetwork.OXYGEN_TENSION: 0.6})
        o2_at_1_0 = self.event.calculate_rate_at_patch(self.network, 1)
        self.assertEqual(o2_at_1_0, self.rp * 2 * ((1.0 ** self.sig) / (1.0 ** self.sig + self.hs ** self.sig)))

        # Increase in oxygen increases the rate of change
        self.assertTrue(o2_at_1_0 > o2_at_0_4)

    def test_perform(self):
        self.network.update_patch(1, {BACTERIUM_EXTRACELLULAR_DORMANT: 1})
        self.event.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_EXTRACELLULAR_REPLICATING), 1)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_EXTRACELLULAR_DORMANT), 0)


if __name__ == '__main__':
    unittest.main()
