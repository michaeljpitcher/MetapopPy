import unittest
from metapoppy_group_by_patch import *

compartments = ['a','b','c']
attributes = ['d','e','f']


class NAEvent(Event):
    REACTION_PARAMETER = 0.5

    def __init__(self):
        Event.__init__(self)

    def _calculate_state_variable(self, patch, edges):
        return patch[Network.COMPARTMENTS][compartments[0]]

    def perform(self, network):
        network.update_patch(self._patch_id, {compartments[1]: 2})


class EventTestCase(unittest.TestCase):

    def setUp(self):
        self.network = Network(compartments, attributes, [])
        self.network.add_node(1)
        self.network.prepare()

        self.event = NAEvent()
        self.network.attach_event(self.event, 1)

    def test_update_rate(self):
        self.assertFalse(self.event.rate())
        self.event.update_rate(self.network.node[1], [])
        self.assertFalse(self.event.rate())
        # Going directly into the network to update patch
        self.network.node[1][Network.COMPARTMENTS][compartments[0]] = 1
        self.event.update_rate(self.network.node[1], self.network.edges([1]))
        self.assertEqual(self.event.rate(), NAEvent.REACTION_PARAMETER * 1)

        self.network.node[1][Network.COMPARTMENTS][compartments[0]] = 4
        self.event.update_rate(self.network.node[1], self.network.edges([1]))
        self.assertEqual(self.event.rate(), NAEvent.REACTION_PARAMETER * 4)


if __name__ == '__main__':
    unittest.main()
