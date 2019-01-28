import unittest
from tbmetapoppy import *


class PulmonaryNetworkTestCase(unittest.TestCase):

    def setUp(self):
        self.boundary = [(0,5), (0,10), (10,10), (10,0), (0,0)]
        self.network_config = {TBPulmonaryNetwork.TOPOLOGY: TBPulmonaryNetwork.SPACE_FILLING_TREE_2D,
                               TBPulmonaryNetwork.BOUNDARY: self.boundary,
                               TBPulmonaryNetwork.LENGTH_DIVISOR: 2,
                               TBPulmonaryNetwork.MINIMUM_AREA: 6}
        self.compartments = ['a', 'b', 'c']
        self.tree_network = TBPulmonaryNetwork(self.network_config)
        self.tree_network.reset()

    def test_2d_space_filling_tree_initialise(self):
        self.assertEqual(len(self.tree_network.nodes()), 17)
        self.assertEqual(len(self.tree_network.get_patches_by_type(TBPulmonaryNetwork.ALVEOLAR_PATCH)), 16)
        self.assertEqual(len(self.tree_network.get_patches_by_type(TBPulmonaryNetwork.LYMPH_PATCH)), 1)

        positions = [1.25, 3.75, 6.25, 8.75]
        e = []
        for q in positions:
            for w in positions:
                e.append((q,w))
        self.assertItemsEqual(e, self.tree_network._alveolar_positions.values())
        self.assertEqual(len(self.tree_network.edges), 16)
        self.assertEqual(len(self.tree_network[TBPulmonaryNetwork.LYMPH_PATCH]), 16)

    def test_pulmonary_attribute_seeding(self):
        params = {TBPulmonaryNetwork.VENTILATION_SKEW: 3.5, TBPulmonaryNetwork.PERFUSION_SKEW: 3,
                  TBPulmonaryNetwork.DRAINAGE_SKEW: 2}
        seeding = self.tree_network.get_pulmonary_att_seeding(params)

        alv_patch_ids = self.tree_network.get_patches_by_type(TBPulmonaryNetwork.ALVEOLAR_PATCH)

        for att in [TBPulmonaryNetwork.VENTILATION, TBPulmonaryNetwork.PERFUSION]:
            self.assertEqual(1, sum(seeding[a][TypedNetwork.ATTRIBUTES][att] for a in alv_patch_ids))

        maxv = max([seeding[a][TypedNetwork.ATTRIBUTES][TBPulmonaryNetwork.VENTILATION] for a in alv_patch_ids])
        minv = min([seeding[a][TypedNetwork.ATTRIBUTES][TBPulmonaryNetwork.VENTILATION] for a in alv_patch_ids])

        self.assertAlmostEqual(maxv / minv, params[TBPulmonaryNetwork.VENTILATION_SKEW])

        maxq = max([seeding[a][TypedNetwork.ATTRIBUTES][TBPulmonaryNetwork.PERFUSION] for a in alv_patch_ids])
        minq = min([seeding[a][TypedNetwork.ATTRIBUTES][TBPulmonaryNetwork.PERFUSION] for a in alv_patch_ids])

        self.assertAlmostEqual(maxq / minq, params[TBPulmonaryNetwork.PERFUSION_SKEW])

        maxd = max([seeding[a][TypedNetwork.ATTRIBUTES][TBPulmonaryNetwork.DRAINAGE] for a in alv_patch_ids])
        mind = min([seeding[a][TypedNetwork.ATTRIBUTES][TBPulmonaryNetwork.DRAINAGE] for a in alv_patch_ids])

        self.assertEqual(mind, 2.0 / 3.0)
        self.assertEqual(maxd, 4.0 / 3.0)

        self.assertAlmostEqual(maxd / mind, params[TBPulmonaryNetwork.DRAINAGE_SKEW])

        for a in alv_patch_ids:
            self.assertEqual(seeding[a][TypedNetwork.ATTRIBUTES][TBPulmonaryNetwork.OXYGEN_TENSION],
                             seeding[a][TypedNetwork.ATTRIBUTES][TBPulmonaryNetwork.VENTILATION] /
                             seeding[a][TypedNetwork.ATTRIBUTES][TBPulmonaryNetwork.PERFUSION])


if __name__ == '__main__':
    unittest.main()
