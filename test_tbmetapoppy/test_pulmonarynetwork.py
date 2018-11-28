import unittest
from tbmetapoppy import *


class PulmonaryNetworkTestCase(unittest.TestCase):

    def setUp(self):
        self.boundary = [(0,5), (0,10), (10,10), (10,0), (0,0)]
        self.network_config = {PulmonaryNetwork.TOPOLOGY: PulmonaryNetwork.SPACE_FILLING_TREE_2D,
                               PulmonaryNetwork.BOUNDARY: self.boundary,
                               PulmonaryNetwork.LENGTH_DIVISOR: 2,
                               PulmonaryNetwork.MINIMUM_AREA: 6}
        self.compartments = ['a', 'b', 'c']
        self.tree_network = PulmonaryNetwork(self.network_config, self.compartments)
        self.tree_network.reset()

    def test_2d_space_filling_tree_initialise(self):
        self.assertEqual(len(self.tree_network.nodes()), 17)
        self.assertEqual(len(self.tree_network.get_patches_by_type(PulmonaryNetwork.ALVEOLAR_PATCH)), 16)
        self.assertEqual(len(self.tree_network.get_patches_by_type(PulmonaryNetwork.LYMPH_PATCH)), 1)

        positions = [1.25, 3.75, 6.25, 8.75]
        e = []
        for q in positions:
            for w in positions:
                e.append((q,w))
        self.assertItemsEqual(e, self.tree_network._alveolar_positions.values())
        self.assertEqual(len(self.tree_network.edges), 16)
        self.assertEqual(len(self.tree_network[PulmonaryNetwork.LYMPH_PATCH]), 16)

    def test_pulmonary_attribute_seeding(self):
        params = {TBDynamics.VENTILATION_SKEW: 3.5, TBDynamics.PERFUSION_SKEW: 3,
                  TBDynamics.DRAINAGE_SKEW: 2}
        seeding = self.tree_network.get_pulmonary_att_seeding(params[TBDynamics.VENTILATION_SKEW],
                                                              params[TBDynamics.PERFUSION_SKEW],
                                                              params[TBDynamics.DRAINAGE_SKEW])

        alv_patch_ids = self.tree_network.get_patches_by_type(PulmonaryNetwork.ALVEOLAR_PATCH)

        for att in [PulmonaryNetwork.VENTILATION, PulmonaryNetwork.PERFUSION]:
            self.assertEqual(1, sum(seeding[a][TypedNetwork.ATTRIBUTES][att] for a in alv_patch_ids))

        maxv = max([seeding[a][TypedNetwork.ATTRIBUTES][PulmonaryNetwork.VENTILATION] for a in alv_patch_ids])
        minv = min([seeding[a][TypedNetwork.ATTRIBUTES][PulmonaryNetwork.VENTILATION] for a in alv_patch_ids])

        self.assertAlmostEqual(maxv / minv, params[TBDynamics.VENTILATION_SKEW])

        maxq = max([seeding[a][TypedNetwork.ATTRIBUTES][PulmonaryNetwork.PERFUSION] for a in alv_patch_ids])
        minq = min([seeding[a][TypedNetwork.ATTRIBUTES][PulmonaryNetwork.PERFUSION] for a in alv_patch_ids])

        self.assertAlmostEqual(maxq / minq, params[TBDynamics.PERFUSION_SKEW])

        maxd = max([seeding[a][TypedNetwork.ATTRIBUTES][PulmonaryNetwork.DRAINAGE] for a in alv_patch_ids])
        mind = min([seeding[a][TypedNetwork.ATTRIBUTES][PulmonaryNetwork.DRAINAGE] for a in alv_patch_ids])

        self.assertEqual(mind, 2.0 / 3.0)
        self.assertEqual(maxd, 4.0 / 3.0)

        self.assertAlmostEqual(maxd / mind, params[TBDynamics.DRAINAGE_SKEW])

        for a in alv_patch_ids:
            self.assertEqual(seeding[a][TypedNetwork.ATTRIBUTES][PulmonaryNetwork.OXYGEN_TENSION],
                             seeding[a][TypedNetwork.ATTRIBUTES][PulmonaryNetwork.VENTILATION] /
                             seeding[a][TypedNetwork.ATTRIBUTES][PulmonaryNetwork.PERFUSION])


if __name__ == '__main__':
    unittest.main()
