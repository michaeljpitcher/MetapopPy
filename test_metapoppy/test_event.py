import unittest
from metapoppy import *

compartments = ['a','b','c']
attributes = ['d','e','f']

RP1_key = 'rp1_key'
RP2_key = 'rp2_key'
PARAM_A = 'param_a'
PARAM_B = 'param_b'
PARAM_C = 'param_c'


class EventTestCase(unittest.TestCase):

    def setUp(self):
        self.network = Network(compartments, attributes, [])
        self.network.add_node(1)
        self.network.prepare()

        class NAEvent(Event):
            def __init__(self, rp_key):
                Event.__init__(self, rp_key)

            def _calculate_state_variable_at_patch(self, network, patch_id):
                return network.get_compartment_value(patch_id, compartments[0])

            def perform(self, network, patch_id):
                network.update_patch(patch_id, {compartments[1]: 1})

        self.event = NAEvent(RP1_key)

    def test_calculate_rate(self):
        self.event.set_parameters({RP1_key: 0.1})
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
        self.typed_attributes = {self.patch_types[0]: attributes[0], self.patch_types[1]: attributes[1]}
        self.network = TypedNetwork(compartments, self.typed_attributes, [])
        self.network.add_node(1)
        self.network.set_patch_type(1, self.patch_types[0])
        self.network.add_node(2)
        self.network.set_patch_type(2, self.patch_types[1])
        self.network.prepare()

        class NAPatchTypeEvent(PatchTypeEvent):
            def __init__(self, patch_type, rp_key, add_par1, add_par2):
                self.add_par1 = add_par1
                self.add_par2 = add_par2
                PatchTypeEvent.__init__(self, patch_type, rp_key, [add_par1, add_par2])

            def _calculate_state_variable_at_patch(self, network, patch_id):
                return network.get_compartment_value(patch_id, compartments[0]) * (self._parameters[self.add_par1] +
                                                                                   self._parameters[self.add_par2])

            def perform(self, network, patch_id):
                network.update_patch(patch_id, {compartments[1]: 1})

        self.event_type1 = NAPatchTypeEvent(self.patch_types[0], RP1_key, PARAM_A, PARAM_B)
        self.event_type2 = NAPatchTypeEvent(self.patch_types[1], RP2_key, PARAM_A, PARAM_C)

    def test_calculate_rate(self):
        params = {RP1_key: 0.1, RP2_key: 0.2, PARAM_A: 2, PARAM_B: 3, PARAM_C:5}
        self.event_type1.set_parameters(params)
        self.event_type2.set_parameters(params)

        self.assertFalse(self.event_type1.calculate_rate_at_patch(self.network, 1))
        self.assertFalse(self.event_type1.calculate_rate_at_patch(self.network, 2))
        self.assertFalse(self.event_type2.calculate_rate_at_patch(self.network, 1))
        self.assertFalse(self.event_type2.calculate_rate_at_patch(self.network, 2))

        self.network.node[1][Network.COMPARTMENTS][compartments[0]] = 10

        self.assertAlmostEqual(self.event_type1.calculate_rate_at_patch(self.network, 1), params[RP1_key] * 10 *
                         (params[PARAM_A] + params[PARAM_B]))
        self.assertFalse(self.event_type2.calculate_rate_at_patch(self.network, 2))

        self.network.node[2][Network.COMPARTMENTS][compartments[0]] = 9

        self.assertFalse(self.event_type1.calculate_rate_at_patch(self.network, 2))
        self.assertAlmostEqual(self.event_type2.calculate_rate_at_patch(self.network, 2), params[RP2_key] * 9 *
                         (params[PARAM_A] + params[PARAM_C]))


if __name__ == '__main__':
    unittest.main()
