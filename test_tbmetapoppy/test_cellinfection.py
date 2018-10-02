import unittest
from tbmetapoppy import *


class CellInfectionTestCase(unittest.TestCase):

    def setUp(self):
        self.event_mac = CellInfection(MACROPHAGE_RESTING)
        self.event_dc = CellInfection(DENDRITIC_CELL_IMMATURE)
        self.rp = 0.1
        self.hs1 = 100
        self.hs2 = 88
        self.event_mac.set_parameters(self.rp, self.hs1)
        self.event_dc.set_parameters(self.rp, self.hs2)
        self.network = PulmonaryNetwork({PulmonaryNetwork.TOPOLOGY: None}, TB_COMPARTMENTS)
        self.network.add_node(1)
        self.network.set_patch_type(1, PulmonaryNetwork.ALVEOLAR_PATCH)
        self.network.prepare()

    def test_rate(self):
        self.assertFalse(self.event_mac.calculate_rate_at_patch(self.network, 1))
        self.assertFalse(self.event_dc.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, compartment_changes={MACROPHAGE_RESTING: 9})
        self.assertFalse(self.event_mac.calculate_rate_at_patch(self.network, 1))
        self.assertFalse(self.event_dc.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, compartment_changes={BACTERIUM_EXTRACELLULAR_REPLICATING: 11})
        self.assertEqual(self.event_mac.calculate_rate_at_patch(self.network, 1), self.rp * 9 * (11.0 / (11 + self.hs1)))
        self.assertFalse(self.event_dc.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, compartment_changes={BACTERIUM_EXTRACELLULAR_DORMANT: 6})
        self.assertEqual(self.event_mac.calculate_rate_at_patch(self.network, 1), self.rp * 9 * (17.0 / (17 + self.hs1)))
        self.assertFalse(self.event_dc.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, compartment_changes={DENDRITIC_CELL_IMMATURE: 3})
        self.assertEqual(self.event_mac.calculate_rate_at_patch(self.network, 1), self.rp * 9 * (17.0 / (17 + self.hs1)))
        self.assertEqual(self.event_dc.calculate_rate_at_patch(self.network, 1), self.rp * 3 * (17.0 / (17 + self.hs2)))

    def test_perform(self):
        self.network.update_patch(1, {MACROPHAGE_RESTING: 101,
                                      BACTERIUM_EXTRACELLULAR_REPLICATING: 100,
                                      BACTERIUM_EXTRACELLULAR_DORMANT: 1})
        for e in range(1,102):
            self.event_mac.perform(self.network, 1)
            self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_EXTRACELLULAR_DORMANT) +
                             self.network.get_compartment_value(1, BACTERIUM_EXTRACELLULAR_REPLICATING), 101-e)

        self.assertEqual(self.network.get_compartment_value(1, MACROPHAGE_RESTING), 0)
        self.assertEqual(self.network.get_compartment_value(1, MACROPHAGE_INFECTED), 101)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_EXTRACELLULAR_DORMANT), 0)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_EXTRACELLULAR_REPLICATING), 0)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_INTRACELLULAR_MACROPHAGE), 101)


if __name__ == '__main__':
    unittest.main()
