import unittest
from tbmetapoppy import *
from metapoppy.environment import Environment


class ReplicationTestCase(unittest.TestCase):

    def setUp(self):
        self.event = Replication(TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_REPLICATING)
        self.rep_key = self.event.parameter_keys()[0]
        self.params = {self.rep_key: 0.1}
        self.event.set_parameters(self.params)
        self.network = TBPulmonaryEnvironment({TBPulmonaryEnvironment.TOPOLOGY: None})
        self.network.add_edge(1, TBPulmonaryEnvironment.LYMPH_PATCH)
        self.network.set_patch_type(1, TBPulmonaryEnvironment.ALVEOLAR_PATCH)
        self.network.set_patch_type(TBPulmonaryEnvironment.LYMPH_PATCH, TBPulmonaryEnvironment.LYMPH_PATCH)
        self.network.reset()

    def test_initialise(self):
        self.assertItemsEqual([TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_REPLICATING + Replication.REPLICATION_RATE],
                              self.event.parameter_keys())

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_REPLICATING: 1})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.params[self.rep_key] * 1)
        self.network.update_patch(1, {TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_REPLICATING: 3})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), self.params[self.rep_key] * 4)

    def test_perform(self):
        self.network.update_patch(1, {TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_REPLICATING: 1})
        self.event.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_REPLICATING), 2)


class IntracellularBacterialReplicationTestCase(unittest.TestCase):

    def setUp(self):
        self.event = IntracellularBacterialReplication()

        self.rep_key = TBPulmonaryEnvironment.BACTERIUM_INTRACELLULAR_MACROPHAGE + Replication.REPLICATION_RATE

        self.params = {self.rep_key: 0.1, INTRACELLULAR_REPLICATION_SIGMOID: 2, MACROPHAGE_CAPACITY: 30}
        self.event.set_parameters(self.params)
        self.network = TBPulmonaryEnvironment({TBPulmonaryEnvironment.TOPOLOGY: None})
        self.network.add_edge(1, TBPulmonaryEnvironment.LYMPH_PATCH)
        self.network.set_patch_type(1, TBPulmonaryEnvironment.ALVEOLAR_PATCH)
        self.network.set_patch_type(TBPulmonaryEnvironment.LYMPH_PATCH, TBPulmonaryEnvironment.LYMPH_PATCH)
        self.network.reset()

    def test_initialise(self):
        self.assertItemsEqual(self.params.keys(), self.event.parameter_keys())

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, {TBPulmonaryEnvironment.BACTERIUM_INTRACELLULAR_MACROPHAGE: 18,
                                      TBPulmonaryEnvironment.MACROPHAGE_INFECTED: 5})



        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), 0.1 * 18 * (1 - (
                    1.0 * (18 ** self.params[INTRACELLULAR_REPLICATION_SIGMOID]) / (
                        18 ** self.params[INTRACELLULAR_REPLICATION_SIGMOID] + (5 * self.params[MACROPHAGE_CAPACITY])
                        ** self.params[INTRACELLULAR_REPLICATION_SIGMOID]))))

    def test_perform(self):
        self.network.update_patch(1, {TBPulmonaryEnvironment.MACROPHAGE_INFECTED: 1,
                                      TBPulmonaryEnvironment.BACTERIUM_INTRACELLULAR_MACROPHAGE: 10})
        self.event.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryEnvironment.BACTERIUM_INTRACELLULAR_MACROPHAGE), 11)


if __name__ == '__main__':
    unittest.main()
