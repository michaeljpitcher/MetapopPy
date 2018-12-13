import unittest
from tbmetapoppy.events.cell_translocation import *

TEST_TRANSLOCATION_RATE_MI = 'test_translocation_rate_mi'
TEST_TRANSLOCATION_RATE_DCM = 'test_translocation_rate_dcm'
TEST_TRANSLOCATION_RATE_TA = 'test_translocation_rate_ta'


class CellTranslocationToLymphTestCase(unittest.TestCase):

    def setUp(self):
        self.event_mi = CellTranslocationToLymph(TEST_TRANSLOCATION_RATE_MI, TBPulmonaryNetwork.MACROPHAGE_INFECTED)
        self.event_dcm = CellTranslocationToLymph(TEST_TRANSLOCATION_RATE_DCM, TBPulmonaryNetwork.DENDRITIC_CELL_MATURE)
        self.params = {TEST_TRANSLOCATION_RATE_MI: 0.1, TEST_TRANSLOCATION_RATE_DCM: 0.2}
        self.event_mi.set_parameters(self.params)
        self.event_dcm.set_parameters(self.params)
        self.network = TBPulmonaryNetwork({TBPulmonaryNetwork.TOPOLOGY: None})
        self.network.add_edge(1, TBPulmonaryNetwork.LYMPH_PATCH)
        self.network.set_patch_type(1, TBPulmonaryNetwork.ALVEOLAR_PATCH)
        self.network.set_patch_type(TBPulmonaryNetwork.LYMPH_PATCH, TBPulmonaryNetwork.LYMPH_PATCH)
        self.network.reset()

    def test_rate(self):
        self.assertFalse(self.event_mi.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, compartment_changes={TBPulmonaryNetwork.MACROPHAGE_INFECTED: 3, TBPulmonaryNetwork.DENDRITIC_CELL_MATURE: 5},
                                  attribute_changes={TBPulmonaryNetwork.DRAINAGE: 0.7})
        self.assertEqual(self.event_mi.calculate_rate_at_patch(self.network, 1),
                         self.params[TEST_TRANSLOCATION_RATE_MI] * 0.7 * 3)
        self.assertAlmostEqual(self.event_dcm.calculate_rate_at_patch(self.network, 1),
                         self.params[TEST_TRANSLOCATION_RATE_DCM] * 0.7 * 5)

    def test_perform(self):
        # MI even split
        self.network.update_patch(1, compartment_changes={TBPulmonaryNetwork.MACROPHAGE_INFECTED: 2,
                                                          TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE: 10})
        self.event_mi.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.MACROPHAGE_INFECTED), 1)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE), 5)
        self.assertEqual(self.network.get_compartment_value(TBPulmonaryNetwork.LYMPH_PATCH,
                                                            TBPulmonaryNetwork.MACROPHAGE_INFECTED), 1)
        self.assertEqual(self.network.get_compartment_value(TBPulmonaryNetwork.LYMPH_PATCH,
                                                            TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE), 5)

        # MI round down
        self.network.update_patch(1, compartment_changes={TBPulmonaryNetwork.MACROPHAGE_INFECTED: 2,
                                                          TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE: 5})
        self.network.update_patch(TBPulmonaryNetwork.LYMPH_PATCH, compartment_changes={TBPulmonaryNetwork.MACROPHAGE_INFECTED: -1,
                                                                                       TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE: -5})
        self.event_mi.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.MACROPHAGE_INFECTED), 2)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE), 7)
        self.assertEqual(self.network.get_compartment_value(TBPulmonaryNetwork.LYMPH_PATCH,
                                                            TBPulmonaryNetwork.MACROPHAGE_INFECTED), 1)
        self.assertEqual(self.network.get_compartment_value(TBPulmonaryNetwork.LYMPH_PATCH,
                                                            TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE), 3)

        # MI round up
        self.network.update_patch(1, compartment_changes={TBPulmonaryNetwork.MACROPHAGE_INFECTED: 2,
                                                          TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE: 3})
        self.network.update_patch(TBPulmonaryNetwork.LYMPH_PATCH,
                                  compartment_changes={TBPulmonaryNetwork.MACROPHAGE_INFECTED: -1,
                                                       TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE: -3})
        self.event_mi.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.MACROPHAGE_INFECTED), 3)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE), 7)
        self.assertEqual(self.network.get_compartment_value(TBPulmonaryNetwork.LYMPH_PATCH,
                                                            TBPulmonaryNetwork.MACROPHAGE_INFECTED), 1)
        self.assertEqual(self.network.get_compartment_value(TBPulmonaryNetwork.LYMPH_PATCH,
                                                            TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE), 3)


class CellTranslocationToLungTestCase(unittest.TestCase):

    def setUp(self):
        self.event = CellTranslocationToLung(TEST_TRANSLOCATION_RATE_TA, TBPulmonaryNetwork.T_CELL_ACTIVATED)
        self.params = {TEST_TRANSLOCATION_RATE_TA: 0.1}
        self.event.set_parameters(self.params)
        self.network = TBPulmonaryNetwork({TBPulmonaryNetwork.TOPOLOGY: None})
        for a in range(1,4):
            self.network.add_edge(a, TBPulmonaryNetwork.LYMPH_PATCH)
            self.network.set_patch_type(a, TBPulmonaryNetwork.ALVEOLAR_PATCH)
        self.network.set_patch_type(TBPulmonaryNetwork.LYMPH_PATCH, TBPulmonaryNetwork.LYMPH_PATCH)
        self.network.reset()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.LYMPH_PATCH))
        self.network.update_patch(TBPulmonaryNetwork.LYMPH_PATCH, {TBPulmonaryNetwork.T_CELL_ACTIVATED: 5})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, TBPulmonaryNetwork.LYMPH_PATCH),
                         self.params[TEST_TRANSLOCATION_RATE_TA] * 5)

    def test_perform(self):
        self.network.update_edge(1, TBPulmonaryNetwork.LYMPH_PATCH, {TBPulmonaryNetwork.PERFUSION: 0.1})
        self.network.update_edge(2, TBPulmonaryNetwork.LYMPH_PATCH, {TBPulmonaryNetwork.PERFUSION: 0.2})
        self.network.update_edge(3, TBPulmonaryNetwork.LYMPH_PATCH, {TBPulmonaryNetwork.PERFUSION: 0.3})

        seed = 10000
        self.network.update_patch(TBPulmonaryNetwork.LYMPH_PATCH, {TBPulmonaryNetwork.T_CELL_ACTIVATED: seed})
        for n in range(seed):
            self.event.perform(self.network, TBPulmonaryNetwork.LYMPH_PATCH)

        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryNetwork.T_CELL_ACTIVATED) + \
              self.network.get_compartment_value(2, TBPulmonaryNetwork.T_CELL_ACTIVATED) + \
              self.network.get_compartment_value(3, TBPulmonaryNetwork.T_CELL_ACTIVATED), seed)

        # TODO - not a true test due to probability of edge choice so could fail
        self.assertTrue(self.network.get_compartment_value(1, TBPulmonaryNetwork.T_CELL_ACTIVATED) <
                        self.network.get_compartment_value(2, TBPulmonaryNetwork.T_CELL_ACTIVATED) <
                        self.network.get_compartment_value(3, TBPulmonaryNetwork.T_CELL_ACTIVATED))



if __name__ == '__main__':
    unittest.main()
