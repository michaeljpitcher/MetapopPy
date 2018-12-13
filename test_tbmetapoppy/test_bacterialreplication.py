import unittest
from tbmetapoppy import *
from metapoppy.network import Network


TEST_REPLICATION_RATE = 'test_replication_rate'
TEST_SIGMOID = 'test_sigmoid'
TEST_HALFSAT = 'test_halfsat'


class ReplicationTestCase(unittest.TestCase):

    def setUp(self):
        self.event = Replication(TEST_REPLICATION_RATE, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING)
        self.params = {TEST_REPLICATION_RATE: 0.1}
        self.event.set_parameters(self.params)
        self.network = TBPulmonaryNetwork({TBPulmonaryNetwork.TOPOLOGY: None})
        self.network.add_edge(1, TBPulmonaryNetwork.LYMPH_PATCH)
        self.network.set_patch_type(1, TBPulmonaryNetwork.ALVEOLAR_PATCH)
        self.network.set_patch_type(TBPulmonaryNetwork.LYMPH_PATCH, TBPulmonaryNetwork.LYMPH_PATCH)
        self.network.reset()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING: 1})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.params[TEST_REPLICATION_RATE] * 1)
        self.network.update_patch(1, {TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING: 3})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.params[TEST_REPLICATION_RATE] * 4)

    def test_perform(self):
        self.network.update_patch(1, {TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING: 1})
        self.event.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING), 2)


class IntracellularBacterialReplicationTestCase(unittest.TestCase):

    def setUp(self):
        self.event = IntracellularBacterialReplication(TEST_REPLICATION_RATE, TEST_SIGMOID, TEST_HALFSAT)
        self.params = {TEST_REPLICATION_RATE: 0.1, TEST_SIGMOID: 2, TEST_HALFSAT: 30}
        self.event.set_parameters(self.params)
        self.network = TBPulmonaryNetwork({TBPulmonaryNetwork.TOPOLOGY: None})
        self.network.add_edge(1, TBPulmonaryNetwork.LYMPH_PATCH)
        self.network.set_patch_type(1, TBPulmonaryNetwork.ALVEOLAR_PATCH)
        self.network.set_patch_type(TBPulmonaryNetwork.LYMPH_PATCH, TBPulmonaryNetwork.LYMPH_PATCH)
        self.network.reset()

    def test_rate(self):

        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE: 18,
                                      TBPulmonaryNetwork.MACROPHAGE_INFECTED: 5})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), 0.1 * 18 * (1 - (
                    1.0 * (18 ** self.params[TEST_SIGMOID]) / (
                        18 ** self.params[TEST_SIGMOID] + (5 * self.params[TEST_HALFSAT]) ** self.params[
                    TEST_SIGMOID]))))

    def test_perform(self):
        self.network.update_patch(1, {TBPulmonaryNetwork.MACROPHAGE_INFECTED: 1,
                                      TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE: 10})
        self.event.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE), 11)


if __name__ == '__main__':
    unittest.main()
