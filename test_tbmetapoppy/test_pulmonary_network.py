import unittest
from tbmetapoppy import *


class PulmonaryNetworkTestCase(unittest.TestCase):

    def setUp(self):
        self.network = PulmonaryNetwork()

    def test_add_alveolar_patch(self):
        self.network.add_alveolar_patch(775)
        self.network.prepare()
        self.assertItemsEqual(self.network.get_patches_by_type(PulmonaryNetwork.ALVEROLAR_PATCH), [775])

    def test_add_lymphatic_patch(self):
        self.network.add_lymphatic_patch(629)
        self.network.prepare()
        self.assertItemsEqual(self.network.get_patches_by_type(PulmonaryNetwork.LYMPHATIC_PATCH), [629])


if __name__ == '__main__':
    unittest.main()
