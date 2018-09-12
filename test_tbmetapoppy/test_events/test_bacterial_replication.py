import unittest
from tbmetapoppy import *


class ExtracellularBacterialReplicationTestCase(unittest.TestCase):

    def setUp(self):
        self.event_ber = ExtracellularBacterialReplication(PulmonaryNetwork.BACTERIUM_EXTRACELLUAR_REPLICATING)
        self.event_bed = ExtracellularBacterialReplication(PulmonaryNetwork.BACTERIUM_EXTRACELLUAR_DORMANT)

    def test_calculate_rate(self):
        patch = {Network.COMPARTMENTS: {PulmonaryNetwork.BACTERIUM_EXTRACELLUAR_REPLICATING: 1,
                                        PulmonaryNetwork.BACTERIUM_EXTRACELLUAR_DORMANT: 2}}



        self.event_ber.calculate_rate(patch, [])
        print self.event_ber.rate()


if __name__ == '__main__':
    unittest.main()
