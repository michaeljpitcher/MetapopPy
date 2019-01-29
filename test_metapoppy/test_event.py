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
        self.network.reset()

        class NAEvent(Event):
            def __init__(self):
                Event.__init__(self, compartments[0], [], [])

            def _define_parameter_keys(self):
                return RP1_key, []

            def _calculate_state_variable_at_patch(self, network, patch_id):
                return network.get_compartment_value(patch_id, compartments[0])

            def perform(self, network, patch_id):
                network.update_patch(patch_id, {compartments[1]: 1})

        self.event = NAEvent()

    def test_calculate_rate(self):
        self.event.set_parameters({RP1_key: 0.1})
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))
        self.network.node[1][Network.COMPARTMENTS][compartments[0]] = 10
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1), 0.1 * 10)

    def test_perform(self):
        self.assertEqual(self.network.get_compartment_value(1, compartments[1]), 0)
        self.event.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, compartments[1]), 1)


class NAPatchTypeEvent(PatchTypeEvent):
    PAR1 = 'par1'
    PAR2 = 'par2'

    def __init__(self, patch_type):
        PatchTypeEvent.__init__(self, patch_type, compartments[0], [], [])

    def _define_parameter_keys(self):
        return "test_" + self._patch_type, [NAPatchTypeEvent.PAR1, NAPatchTypeEvent.PAR2]

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, compartments[0]) * (self._parameters[NAPatchTypeEvent.PAR1] +
                                                                           self._parameters[NAPatchTypeEvent.PAR2])

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {compartments[1]: 1})

class PatchTypeEventTestCase(unittest.TestCase):

    def setUp(self):
        self.patch_types = ['type1', 'type2']
        self.typed_attributes = {self.patch_types[0]: attributes[0], self.patch_types[1]: attributes[1]}
        self.network = TypedNetwork(compartments, self.typed_attributes, [])
        self.network.add_node(1)
        self.network.set_patch_type(1, self.patch_types[0])
        self.network.add_node(2)
        self.network.set_patch_type(2, self.patch_types[1])
        self.network.reset()

        self.event_type1 = NAPatchTypeEvent(self.patch_types[0])
        self.event_type2 = NAPatchTypeEvent(self.patch_types[1])

    def test_calculate_rate(self):
        params = {}
        for p in self.event_type1.parameter_keys():
            params[p] = numpy.random.random()
        for p in self.event_type2.parameter_keys():
            params[p] = numpy.random.random()
        self.event_type1.set_parameters(params)
        self.event_type2.set_parameters(params)

        self.assertFalse(self.event_type1.calculate_rate_at_patch(self.network, 1))
        self.assertFalse(self.event_type1.calculate_rate_at_patch(self.network, 2))
        self.assertFalse(self.event_type2.calculate_rate_at_patch(self.network, 1))
        self.assertFalse(self.event_type2.calculate_rate_at_patch(self.network, 2))

        self.network.node[1][Network.COMPARTMENTS][compartments[0]] = 10

        self.assertAlmostEqual(self.event_type1.calculate_rate_at_patch(self.network, 1), params['test_type1'] * 10 *
                         (params[NAPatchTypeEvent.PAR1] + params[NAPatchTypeEvent.PAR2]))
        self.assertFalse(self.event_type2.calculate_rate_at_patch(self.network, 2))

        self.network.node[2][Network.COMPARTMENTS][compartments[0]] = 9

        self.assertFalse(self.event_type1.calculate_rate_at_patch(self.network, 2))
        self.assertAlmostEqual(self.event_type2.calculate_rate_at_patch(self.network, 2), params['test_type2'] * 9 *
                               (params[NAPatchTypeEvent.PAR1] + params[NAPatchTypeEvent.PAR2]))


if __name__ == '__main__':
    unittest.main()
