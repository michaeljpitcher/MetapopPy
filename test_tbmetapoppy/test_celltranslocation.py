import unittest
from tbmetapoppy.events.cell_translocation import *


class CellTranslocationToLymphTestCase(unittest.TestCase):

    def setUp(self):
        self.event_mi = CellTranslocationToLymph(MACROPHAGE_INFECTED)
        self.event_dcm = CellTranslocationToLymph(DENDRITIC_CELL_MATURE)
        self.rp = 0.1
        self.event_mi.set_reaction_parameter(self.rp)
        self.event_dcm.set_reaction_parameter(self.rp)
        self.network = PulmonaryNetwork({PulmonaryNetwork.TOPOLOGY: None}, TB_COMPARTMENTS)
        self.network.add_edge(1, PulmonaryNetwork.LYMPH_PATCH)
        self.network.set_patch_type(1, PulmonaryNetwork.ALVEOLAR_PATCH)
        self.network.set_patch_type(PulmonaryNetwork.LYMPH_PATCH, PulmonaryNetwork.LYMPH_PATCH)
        self.network.prepare()

    def test_rate(self):
        self.assertFalse(self.event_mi.calculate_rate_at_patch(self.network, 1))
        self.network.update_patch(1, compartment_changes={MACROPHAGE_INFECTED: 3},
                                  attribute_changes={PulmonaryNetwork.DRAINAGE: 0.7})
        self.assertEqual(self.event_mi.calculate_rate_at_patch(self.network, 1), self.rp * 0.7 * 3)

    def test_perform(self):
        # MI even split
        self.network.update_patch(1, compartment_changes={MACROPHAGE_INFECTED: 2,
                                                          BACTERIUM_INTRACELLULAR_MACROPHAGE: 10})
        self.event_mi.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, MACROPHAGE_INFECTED), 1)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_INTRACELLULAR_MACROPHAGE), 5)
        self.assertEqual(self.network.get_compartment_value(PulmonaryNetwork.LYMPH_PATCH,
                                                            MACROPHAGE_INFECTED), 1)
        self.assertEqual(self.network.get_compartment_value(PulmonaryNetwork.LYMPH_PATCH,
                                                            BACTERIUM_INTRACELLULAR_MACROPHAGE), 5)

        # MI round down
        self.network.update_patch(1, compartment_changes={MACROPHAGE_INFECTED: 2,
                                                          BACTERIUM_INTRACELLULAR_MACROPHAGE: 5})
        self.network.update_patch(PulmonaryNetwork.LYMPH_PATCH, compartment_changes={MACROPHAGE_INFECTED: -1,
                                                          BACTERIUM_INTRACELLULAR_MACROPHAGE: -5})
        self.event_mi.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, MACROPHAGE_INFECTED), 2)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_INTRACELLULAR_MACROPHAGE), 7)
        self.assertEqual(self.network.get_compartment_value(PulmonaryNetwork.LYMPH_PATCH,
                                                            MACROPHAGE_INFECTED), 1)
        self.assertEqual(self.network.get_compartment_value(PulmonaryNetwork.LYMPH_PATCH,
                                                            BACTERIUM_INTRACELLULAR_MACROPHAGE), 3)

        # MI round up
        self.network.update_patch(1, compartment_changes={MACROPHAGE_INFECTED: 2,
                                                          BACTERIUM_INTRACELLULAR_MACROPHAGE: 3})
        self.network.update_patch(PulmonaryNetwork.LYMPH_PATCH,
                                  compartment_changes={MACROPHAGE_INFECTED: -1,
                                                       BACTERIUM_INTRACELLULAR_MACROPHAGE: -3})
        self.event_mi.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, MACROPHAGE_INFECTED), 3)
        self.assertEqual(self.network.get_compartment_value(1, BACTERIUM_INTRACELLULAR_MACROPHAGE), 7)
        self.assertEqual(self.network.get_compartment_value(PulmonaryNetwork.LYMPH_PATCH,
                                                            MACROPHAGE_INFECTED), 1)
        self.assertEqual(self.network.get_compartment_value(PulmonaryNetwork.LYMPH_PATCH,
                                                            BACTERIUM_INTRACELLULAR_MACROPHAGE), 3)


class CellTranslocationToLungTestCase(unittest.TestCase):

    def setUp(self):
        self.event = CellTranslocationToLung(T_CELL_ACTIVATED)
        self.rp = 0.1
        self.event.set_reaction_parameter(self.rp)
        self.network = PulmonaryNetwork({PulmonaryNetwork.TOPOLOGY: None}, TB_COMPARTMENTS)
        for a in range(1,4):
            self.network.add_edge(a, PulmonaryNetwork.LYMPH_PATCH)
            self.network.set_patch_type(a, PulmonaryNetwork.ALVEOLAR_PATCH)
        self.network.set_patch_type(PulmonaryNetwork.LYMPH_PATCH, PulmonaryNetwork.LYMPH_PATCH)
        self.network.prepare()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, PulmonaryNetwork.LYMPH_PATCH))
        self.network.update_patch(PulmonaryNetwork.LYMPH_PATCH, {T_CELL_ACTIVATED: 5})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, PulmonaryNetwork.LYMPH_PATCH), self.rp * 5)

    def test_perform(self):
        self.network.update_edge(1, PulmonaryNetwork.LYMPH_PATCH, {PulmonaryNetwork.PERFUSION: 0.1})
        self.network.update_edge(2, PulmonaryNetwork.LYMPH_PATCH, {PulmonaryNetwork.PERFUSION: 0.2})
        self.network.update_edge(3, PulmonaryNetwork.LYMPH_PATCH, {PulmonaryNetwork.PERFUSION: 0.3})

        seed = 10000
        self.network.update_patch(PulmonaryNetwork.LYMPH_PATCH, {T_CELL_ACTIVATED: seed})
        for n in range(seed):
            self.event.perform(self.network, PulmonaryNetwork.LYMPH_PATCH)

        self.assertEqual(self.network.get_compartment_value(1, T_CELL_ACTIVATED) + \
              self.network.get_compartment_value(2, T_CELL_ACTIVATED) + \
              self.network.get_compartment_value(3, T_CELL_ACTIVATED), seed)

        # TODO - not a true test due to probability of edge choice so could fail
        self.assertTrue(self.network.get_compartment_value(1, T_CELL_ACTIVATED) <
                        self.network.get_compartment_value(2, T_CELL_ACTIVATED) <
                        self.network.get_compartment_value(3, T_CELL_ACTIVATED))



if __name__ == '__main__':
    unittest.main()
