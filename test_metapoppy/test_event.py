import unittest
from metapoppy import *

compartments = ['a','b','c']
attributes = ['d','e','f']


class EventTestCase(unittest.TestCase):

    def setUp(self):
        self.network = Network(compartments, attributes, [])
        self.network.add_node(1)
        self.network.prepare()

        class NAEvent(Event):
            def __init__(self):
                Event.__init__(self)

            def _calculate_state_variable_at_patch(self, network, patch_id):
                return network.get_compartment_value(patch_id, compartments[0])

            def perform(self, network, patch_id):
                network.update_patch(patch_id, {compartments[1]:1})

        self.event = NAEvent()

    def test_calculate_rate(self):
        self.event.set_reaction_parameter(0.1)
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.node[1][Network.COMPARTMENTS][compartments[0]] = 10
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), 0.1 * 10)

    def test_perform(self):
        self.assertEqual(self.network.get_compartment_value(1, compartments[1]), 0)
        self.event.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, compartments[1]), 1)


class PatchTypeEventTestCase(unittest.TestCase):

    def setUp(self):
        self.patch_types = ['type1', 'type2']
        self.network = TypedNetwork(compartments, attributes, [])
        self.network.add_node(1)
        self.network.set_patch_type(1, self.patch_types[0])
        self.network.add_node(2)
        self.network.set_patch_type(2, self.patch_types[1])
        self.network.prepare()

        class NAEvent(PatchTypeEvent):
            def __init__(self, patch_type):
                PatchTypeEvent.__init__(self, patch_type)

            def _calculate_state_variable_at_patch(self, network, patch_id):
                return network.get_compartment_value(patch_id, compartments[0])

            def perform(self, network, patch_id):
                network.update_patch(patch_id, {compartments[1]:1})

        self.event_type1 = NAEvent(self.patch_types[0])
        self.event_type2 = NAEvent(self.patch_types[1])

    def test_calculate_rate(self):
        self.event_type1.set_reaction_parameter(0.1)
        self.event_type2.set_reaction_parameter(0.2)

        self.assertFalse(self.event_type1.calculate_rate_at_patch(self.network, 1))
        self.assertFalse(self.event_type1.calculate_rate_at_patch(self.network, 2))
        self.assertFalse(self.event_type2.calculate_rate_at_patch(self.network, 1))
        self.assertFalse(self.event_type2.calculate_rate_at_patch(self.network, 2))

        self.network.node[1][Network.COMPARTMENTS][compartments[0]] = 10

        self.assertEqual(self.event_type1.calculate_rate_at_patch(self.network, 1), 0.1 * 10)
        self.assertFalse(self.event_type2.calculate_rate_at_patch(self.network, 2))

        self.network.node[2][Network.COMPARTMENTS][compartments[0]] = 9

        self.assertFalse(self.event_type1.calculate_rate_at_patch(self.network, 2))
        self.assertEqual(self.event_type2.calculate_rate_at_patch(self.network, 2), 0.2 * 9)


if __name__ == '__main__':
    unittest.main()
