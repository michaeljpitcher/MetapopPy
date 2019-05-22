import unittest
from tbmetapoppy import *


class PulmonaryNetworkTestCase(unittest.TestCase):

    def setUp(self):
        self.boundary = [(0,5), (0,10), (10,10), (10,0), (0,0)]
        self.network_config = {TBPulmonaryEnvironment.TOPOLOGY: TBPulmonaryEnvironment.SPACE_FILLING_TREE_2D,
                               TBPulmonaryEnvironment.BOUNDARY: self.boundary,
                               TBPulmonaryEnvironment.LENGTH_DIVISOR: 2,
                               TBPulmonaryEnvironment.MINIMUM_AREA: 6}
        self.compartments = ['a', 'b', 'c']
        self.tree_network = TBPulmonaryEnvironment(self.network_config)
        self.tree_network.reset()

    def test_2d_space_filling_tree_initialise(self):
        self.assertEqual(len(self.tree_network.nodes()), 17)
        self.assertEqual(len(self.tree_network.get_patches_by_type(TBPulmonaryEnvironment.ALVEOLAR_PATCH)), 16)
        self.assertEqual(len(self.tree_network.get_patches_by_type(TBPulmonaryEnvironment.LYMPH_PATCH)), 1)

        positions = [1.25, 3.75, 6.25, 8.75]
        e = []
        for q in positions:
            for w in positions:
                e.append((q,w))
        self.assertItemsEqual(e, self.tree_network._alveolar_positions.values())
        self.assertEqual(len(self.tree_network.edges), 16)
        self.assertEqual(len(self.tree_network[TBPulmonaryEnvironment.LYMPH_PATCH]), 16)

    def test_pulmonary_attribute_seeding(self):
        params = {TBPulmonaryEnvironment.VENTILATION_SKEW: 3.5, TBPulmonaryEnvironment.PERFUSION_SKEW: 3,
                  TBPulmonaryEnvironment.DRAINAGE_SKEW: 2}
        self.tree_network.calculate_pulmonary_attribute_values(params)
        seeding = self.tree_network._pulmonary_att_seeding

        alv_patch_ids = self.tree_network.get_patches_by_type(TBPulmonaryEnvironment.ALVEOLAR_PATCH)

        for att in [TBPulmonaryEnvironment.VENTILATION, TBPulmonaryEnvironment.PERFUSION]:
            self.assertEqual(1, sum(seeding[a][TypedEnvironment.ATTRIBUTES][att] for a in alv_patch_ids))

        maxv = max([seeding[a][TypedEnvironment.ATTRIBUTES][TBPulmonaryEnvironment.VENTILATION] for a in alv_patch_ids])
        minv = min([seeding[a][TypedEnvironment.ATTRIBUTES][TBPulmonaryEnvironment.VENTILATION] for a in alv_patch_ids])

        self.assertAlmostEqual(maxv / minv, params[TBPulmonaryEnvironment.VENTILATION_SKEW])

        maxq = max([seeding[a][TypedEnvironment.ATTRIBUTES][TBPulmonaryEnvironment.PERFUSION] for a in alv_patch_ids])
        minq = min([seeding[a][TypedEnvironment.ATTRIBUTES][TBPulmonaryEnvironment.PERFUSION] for a in alv_patch_ids])

        self.assertAlmostEqual(maxq / minq, params[TBPulmonaryEnvironment.PERFUSION_SKEW])

        maxd = max([seeding[a][TypedEnvironment.ATTRIBUTES][TBPulmonaryEnvironment.DRAINAGE] for a in alv_patch_ids])
        mind = min([seeding[a][TypedEnvironment.ATTRIBUTES][TBPulmonaryEnvironment.DRAINAGE] for a in alv_patch_ids])

        self.assertEqual(mind, 2.0 / 3.0)
        self.assertEqual(maxd, 4.0 / 3.0)

        self.assertAlmostEqual(maxd / mind, params[TBPulmonaryEnvironment.DRAINAGE_SKEW])

        for a in alv_patch_ids:
            self.assertEqual(seeding[a][TypedEnvironment.ATTRIBUTES][TBPulmonaryEnvironment.OXYGEN_TENSION],
                             seeding[a][TypedEnvironment.ATTRIBUTES][TBPulmonaryEnvironment.VENTILATION] /
                             seeding[a][TypedEnvironment.ATTRIBUTES][TBPulmonaryEnvironment.PERFUSION])


if __name__ == '__main__':
    unittest.main()
