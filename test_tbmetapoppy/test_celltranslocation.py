import unittest
from tbmetapoppy import *


class TranslocationTestCase(unittest.TestCase):

    def setUp(self):
        self.event = Translocation(TBPulmonaryEnvironment.ALVEOLAR_PATCH, TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE)
        self.params = {'d_m_translocation_from_alveolar_patch_rate': 0.1}
        self.event.set_parameters(self.params)

        self.network = TBPulmonaryEnvironment({TBPulmonaryEnvironment.TOPOLOGY: None})
        self.network.add_edge(1, TBPulmonaryEnvironment.LYMPH_PATCH)
        self.network.set_patch_type(1, TBPulmonaryEnvironment.ALVEOLAR_PATCH)
        self.network.set_patch_type(TBPulmonaryEnvironment.LYMPH_PATCH, TBPulmonaryEnvironment.LYMPH_PATCH)
        self.network.reset()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))

        self.network.update_patch(1, {TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE: 2})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1),
                         self.params[self.event.reaction_parameter()] * 2)

        self.network.update_patch(1, {TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE: 13})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1),
                         self.params[self.event.reaction_parameter()] * 15)

    def test_perform(self):
        self.network.update_patch(1, {TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE: 3,
                                      TBPulmonaryEnvironment.BACTERIUM_INTRACELLULAR_DENDRITIC: 3})
        self.event.perform(self.network, 1)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE), 2)
        self.assertEqual(self.network.get_compartment_value(1, TBPulmonaryEnvironment.BACTERIUM_INTRACELLULAR_DENDRITIC), 2)
        self.assertEqual(self.network.get_compartment_value(TBPulmonaryEnvironment.LYMPH_PATCH,
                                                            TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE), 1)
        self.assertEqual(self.network.get_compartment_value(TBPulmonaryEnvironment.LYMPH_PATCH,
                                                            TBPulmonaryEnvironment.BACTERIUM_INTRACELLULAR_DENDRITIC), 1)


class TranslocationLungToLymphTestCase(unittest.TestCase):
    def setUp(self):
        self.event = TranslocationLungToLymph(TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE)
        self.params = {'d_m_translocation_from_alveolar_patch_rate': 0.1}
        self.event.set_parameters(self.params)

        self.network = TBPulmonaryEnvironment({TBPulmonaryEnvironment.TOPOLOGY: None})
        self.network.add_edge(1, TBPulmonaryEnvironment.LYMPH_PATCH)
        self.network.set_patch_type(1, TBPulmonaryEnvironment.ALVEOLAR_PATCH)
        self.network.set_patch_type(TBPulmonaryEnvironment.LYMPH_PATCH, TBPulmonaryEnvironment.LYMPH_PATCH)
        self.network.reset()

    def test_rate(self):
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))

        self.network.update_patch(1, {TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE: 2})
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, 1))

        self.network.update_patch(1, attribute_changes={TBPulmonaryEnvironment.DRAINAGE: 0.9})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1),
                         self.params[self.event.reaction_parameter()] * 2 * 0.9)

        self.network.update_patch(1, {TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE: 13})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, 1),
                         self.params[self.event.reaction_parameter()] * 15 * 0.9)


class TranslocationLymphToLungTestCase(unittest.TestCase):
    def setUp(self):
        self.event = TranslocationLymphToLung(TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT)
        self.params = {'b_ed_translocation_from_lymph_patch_rate': 0.1}
        self.event.set_parameters(self.params)

        self.network = TBPulmonaryEnvironment({TBPulmonaryEnvironment.TOPOLOGY: None})
        for n in range(5):
            self.network.add_edge(n, TBPulmonaryEnvironment.LYMPH_PATCH)
            self.network.set_patch_type(n, TBPulmonaryEnvironment.ALVEOLAR_PATCH)
        self.network.set_patch_type(TBPulmonaryEnvironment.LYMPH_PATCH, TBPulmonaryEnvironment.LYMPH_PATCH)
        self.network.reset()

    def test_rate(self):
        # No bac
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, TBPulmonaryEnvironment.LYMPH_PATCH))

        # Bac
        self.network.update_patch(TBPulmonaryEnvironment.LYMPH_PATCH,
                                  {TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT: 13})
        self.assertEqual(self.event.calculate_rate_at_patch(self.network, TBPulmonaryEnvironment.LYMPH_PATCH),
                         self.params['b_ed_translocation_from_lymph_patch_rate'] * 13)


    def test_perform(self):
        for n in range(5):
            self.network.update_patch(n, attribute_changes={TBPulmonaryEnvironment.PERFUSION: n + 0.5})

        total = 10000
        self.network.update_patch(TBPulmonaryEnvironment.LYMPH_PATCH,
                                  {TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT: total})

        for n in range(total):
            self.event.perform(self.network, TBPulmonaryEnvironment.LYMPH_PATCH)
        self.assertEqual(sum([self.network.get_compartment_value(n, TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT) for n in range(5)]), total)
        # TODO - random involved so may fail
        self.assertTrue(self.network.get_compartment_value(0, TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT) <
                        self.network.get_compartment_value(1, TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT) <
                        self.network.get_compartment_value(2, TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT) <
                        self.network.get_compartment_value(3, TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT) <
                        self.network.get_compartment_value(4, TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT))


class TCellTranslocationLymphToLungTestCase(unittest.TestCase):
    def setUp(self):
        self.event = TCellTranslocationLymphToLung(TBPulmonaryEnvironment.T_CELL_ACTIVATED)
        self.params = {'t_a_translocation_from_lymph_patch_rate': 0.1,
                       't_a_translocation_from_lymph_patch_sigmoid': 2,
                       't_a_translocation_from_lymph_patch_half_sat': 10}
        self.event.set_parameters(self.params)

        self.network = TBPulmonaryEnvironment({TBPulmonaryEnvironment.TOPOLOGY: None})
        for n in range(5):
            self.network.add_edge(n, TBPulmonaryEnvironment.LYMPH_PATCH)
            self.network.set_patch_type(n, TBPulmonaryEnvironment.ALVEOLAR_PATCH)
        self.network.set_patch_type(TBPulmonaryEnvironment.LYMPH_PATCH, TBPulmonaryEnvironment.LYMPH_PATCH)
        self.network.reset()

    def test_rate(self):
        # No t-cells
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, TBPulmonaryEnvironment.LYMPH_PATCH))

        # No t-cells
        self.network.update_patch(1, {TBPulmonaryEnvironment.MACROPHAGE_INFECTED: 1})
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, TBPulmonaryEnvironment.LYMPH_PATCH))

        # T-cells, no DC
        self.network.update_patch(TBPulmonaryEnvironment.LYMPH_PATCH, {TBPulmonaryEnvironment.T_CELL_ACTIVATED: 7})
        self.assertFalse(self.event.calculate_rate_at_patch(self.network, TBPulmonaryEnvironment.LYMPH_PATCH))

        # T-cells, DC
        self.network.update_patch(TBPulmonaryEnvironment.LYMPH_PATCH, {TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE: 13})
        self.assertAlmostEqual(self.event.calculate_rate_at_patch(self.network, TBPulmonaryEnvironment.LYMPH_PATCH),
                               self.params['t_a_translocation_from_lymph_patch_rate'] * 7.0 *
                               ((13.0)**self.params['t_a_translocation_from_lymph_patch_sigmoid']) /
                               ((13.0)**self.params['t_a_translocation_from_lymph_patch_sigmoid'] +
                          self.params['t_a_translocation_from_lymph_patch_half_sat']**
                          self.params['t_a_translocation_from_lymph_patch_sigmoid']))

        self.network.update_patch(TBPulmonaryEnvironment.LYMPH_PATCH, {TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE: 6})
        self.assertAlmostEqual(self.event.calculate_rate_at_patch(self.network, TBPulmonaryEnvironment.LYMPH_PATCH),
                               self.params['t_a_translocation_from_lymph_patch_rate'] * 7.0 *
                               ((19.0) ** self.params['t_a_translocation_from_lymph_patch_sigmoid']) /
                               ((19.0) ** self.params['t_a_translocation_from_lymph_patch_sigmoid'] +
                           self.params['t_a_translocation_from_lymph_patch_half_sat']**
                           self.params['t_a_translocation_from_lymph_patch_sigmoid']))

    def test_perform(self):
        total = 1000
        self.network.update_patch(TBPulmonaryEnvironment.LYMPH_PATCH, {TBPulmonaryEnvironment.T_CELL_ACTIVATED: total})
        self.network.update_patch(0, attribute_changes={TBPulmonaryEnvironment.PERFUSION: 1.0})
        self.network.update_patch(1, compartment_changes={TBPulmonaryEnvironment.MACROPHAGE_INFECTED: 10,
                                                          TBPulmonaryEnvironment.BACTERIUM_INTRACELLULAR_MACROPHAGE: 10},
                                  attribute_changes={TBPulmonaryEnvironment.PERFUSION: 1.0})  # = 10
        self.network.update_patch(2, compartment_changes={TBPulmonaryEnvironment.MACROPHAGE_INFECTED: 40,
                                                          TBPulmonaryEnvironment.BACTERIUM_INTRACELLULAR_MACROPHAGE: 40},
                                  attribute_changes={TBPulmonaryEnvironment.PERFUSION: 0.5})  # = 20
        self.network.update_patch(3, compartment_changes={TBPulmonaryEnvironment.MACROPHAGE_INFECTED: 120,
                                                          TBPulmonaryEnvironment.BACTERIUM_INTRACELLULAR_MACROPHAGE: 120},
                                  attribute_changes={TBPulmonaryEnvironment.PERFUSION: 0.25})  # = 30
        self.network.update_patch(4, compartment_changes={TBPulmonaryEnvironment.MACROPHAGE_INFECTED: 400,
                                                          TBPulmonaryEnvironment.BACTERIUM_INTRACELLULAR_MACROPHAGE: 400},
                                  attribute_changes={TBPulmonaryEnvironment.PERFUSION: 0.1})  # = 40

        for r in range(total):
            self.event.perform(self.network, TBPulmonaryEnvironment.LYMPH_PATCH)

        # Test they all move - stochastic so cannot test exactly where
        self.assertEqual(sum([self.network.get_compartment_value(n, TBPulmonaryEnvironment.T_CELL_ACTIVATED) for n in
                              range(5)]), total)
        self.assertFalse(self.network.get_compartment_value(0, TBPulmonaryEnvironment.T_CELL_ACTIVATED))

        # TODO - not accurate tests as these are probabilities. Have given margin for error but might fail
        self.assertTrue(total * 0.15 > self.network.get_compartment_value(1, TBPulmonaryEnvironment.T_CELL_ACTIVATED) >
                        total * 0.05)
        self.assertTrue(total * 0.25 > self.network.get_compartment_value(2, TBPulmonaryEnvironment.T_CELL_ACTIVATED) >
                        total * 0.15)
        self.assertTrue(total * 0.35 > self.network.get_compartment_value(3, TBPulmonaryEnvironment.T_CELL_ACTIVATED) >
                        total * 0.25)
        self.assertTrue(total * 0.45 > self.network.get_compartment_value(4, TBPulmonaryEnvironment.T_CELL_ACTIVATED) >
                        total * 0.35)
        
if __name__ == '__main__':
    unittest.main()
