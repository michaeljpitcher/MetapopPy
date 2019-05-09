import unittest
from tbmetapoppy import *
import ConfigParser


class TBDynamicsTestCase(unittest.TestCase):

    def setUp(self):
        self.boundary = [(0, 5), (0, 10), (10, 10), (10, 0), (0, 0)]
        network_config = {TBPulmonaryNetwork.TOPOLOGY: TBPulmonaryNetwork.SPACE_FILLING_TREE_2D,
                          TBPulmonaryNetwork.BOUNDARY: self.boundary,
                          TBPulmonaryNetwork.LENGTH_DIVISOR: 2,
                          TBPulmonaryNetwork.MINIMUM_AREA: 6}
        self.dynamics = TBDynamics(network_config)

    def test_initialise(self):
        # Check event creation

        # Bacterial replication
        ber_rep = [e for e in self.dynamics._events if isinstance(e, Replication) and
                   e._cell_type == TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING]
        self.assertEqual(len(ber_rep), 1)
        bed_rep = [e for e in self.dynamics._events if isinstance(e, Replication) and
                   e._cell_type == TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT]
        self.assertEqual(len(bed_rep), 1)
        bim_rep = [e for e in self.dynamics._events if isinstance(e, IntracellularBacterialReplication)]
        self.assertEqual(len(bim_rep), 1)

        # Bacterial state change
        bed_to_ber = [e for e in self.dynamics._events if isinstance(e, BacteriumChangeStateThroughOxygen) and
                      e._compartment_from == TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT]
        self.assertEqual(len(bed_to_ber), 1)
        ber_to_bed = [e for e in self.dynamics._events if isinstance(e, BacteriumChangeStateThroughOxygen) and
                      e._compartment_from == TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING]
        self.assertEqual(len(ber_to_bed), 1)

        # Bacterial translocation
        # bed_transl = [e for e in self.dynamics._events if isinstance(e, CellTranslocationToLung) and
        #               e._cell_type == TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT]
        # self.assertEqual(len(bed_transl), 1)

        # DC recruitment
        dc_rec_lung = [e for e in self.dynamics._events if isinstance(e, StandardCellRecruitmentLung) and
                       e._cell_type == TBPulmonaryNetwork.DENDRITIC_CELL_IMMATURE]
        self.assertEqual(len(dc_rec_lung), 1)
        dc_rec_lung = [e for e in self.dynamics._events if isinstance(e, EnhancedCellRecruitmentLung) and
                       e._cell_type == TBPulmonaryNetwork.DENDRITIC_CELL_IMMATURE]
        self.assertEqual(len(dc_rec_lung), 1)

        # DC infection
        dc_infect = [e for e in self.dynamics._events if isinstance(e, CellIngestBacterium) and
                     e._cell_type == TBPulmonaryNetwork.DENDRITIC_CELL_IMMATURE]
        self.assertEqual(len(dc_infect), 1)

        # DC translocation
        dc_transl = [e for e in self.dynamics._events if isinstance(e, TranslocationLungToLymph) and
                     e._moving_compartment == TBPulmonaryNetwork.DENDRITIC_CELL_MATURE]
        self.assertEqual(len(dc_transl), 1)

        for d in TBPulmonaryNetwork.DENDRITIC_CELLS:
            q = [e for e in self.dynamics._events if isinstance(e, CellDeath) and e._dying_compartment == d]
            self.assertEqual(len(q), 1)

        # Macrophage recruitment
        mac_rec_lung = [e for e in self.dynamics._events if isinstance(e, StandardCellRecruitmentLung) and
                        e._cell_type == TBPulmonaryNetwork.MACROPHAGE_RESTING]
        self.assertEqual(len(mac_rec_lung), 1)
        mac_rec_lymph = [e for e in self.dynamics._events if isinstance(e, StandardCellRecruitmentLymph) and
                         e._cell_type == TBPulmonaryNetwork.MACROPHAGE_RESTING]
        self.assertEqual(len(mac_rec_lymph), 1)

        # Macrophage activation
        mac_act_tc = [e for e in self.dynamics._events if isinstance(e, CellActivation) and
                      e._resting_cell == TBPulmonaryNetwork.MACROPHAGE_RESTING and e._triggers == [TBPulmonaryNetwork.T_CELL_ACTIVATED]]
        self.assertEqual(len(mac_act_tc), 1)
        mac_act_bac = [e for e in self.dynamics._events if isinstance(e, CellActivation) and
                       e._resting_cell == TBPulmonaryNetwork.MACROPHAGE_RESTING and
                       (e._triggers == [TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING]
                        or e._triggers == [TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT])]
        self.assertEqual(len(mac_act_bac), 1)

        # Macrophage death
        for m in [TBPulmonaryNetwork.MACROPHAGE_RESTING, TBPulmonaryNetwork.MACROPHAGE_ACTIVATED]:
            d = [e for e in self.dynamics._events if isinstance(e, CellDeath) and e._dying_compartment == m]
            self.assertEqual(len(d), 1)
        # Not isinstance since bursting and t-cell death are instances (sub-classes)
        mi_nat_death = [e for e in self.dynamics._events if e.__class__ == InfectedCellDeath and
                        e._dying_compartment == TBPulmonaryNetwork.MACROPHAGE_INFECTED]
        self.assertEqual(len(mi_nat_death), 1)
        mi_burst = [e for e in self.dynamics._events if isinstance(e, MacrophageBursting)]
        self.assertEqual(len(mi_burst), 1)
        mi_death_tcell = [e for e in self.dynamics._events if isinstance(e, TCellDestroysMacrophage)]
        self.assertEqual(len(mi_death_tcell), 1)

        # Macrophage ingest bacterium
        mr_ingest_bac = [e for e in self.dynamics._events if isinstance(e, CellIngestBacterium) and
                      e._cell_type == TBPulmonaryNetwork.MACROPHAGE_RESTING]
        self.assertEqual(len(mr_ingest_bac), 1)
        ma_ingest_bac = [e for e in self.dynamics._events if isinstance(e, CellIngestBacterium) and
                         e._cell_type == TBPulmonaryNetwork.MACROPHAGE_ACTIVATED]
        self.assertEqual(len(ma_ingest_bac), 1)

        # Macrophage translocation
        mac_transl = [e for e in self.dynamics._events if isinstance(e, TranslocationLungToLymph) and
                      e._moving_compartment == TBPulmonaryNetwork.MACROPHAGE_INFECTED]
        self.assertEqual(len(mac_transl), 1)

        #  T-cell recruitment
        tc_recruit = [e for e in self.dynamics._events if isinstance(e, StandardCellRecruitmentLymph) and
                      e._cell_type == TBPulmonaryNetwork.T_CELL_NAIVE]
        self.assertEqual(len(tc_recruit), 1)

        # T-cell activation
        tc_act = [e for e in self.dynamics._events if isinstance(e, CellActivation) and
                  e._resting_cell == TBPulmonaryNetwork.T_CELL_NAIVE]
        self.assertEqual(len(tc_act), 1)
        self.assertItemsEqual(tc_act[0]._triggers, [TBPulmonaryNetwork.DENDRITIC_CELL_MATURE, TBPulmonaryNetwork.MACROPHAGE_INFECTED])

        # T-cell translocation
        tc_transl = [e for e in self.dynamics._events if isinstance(e, TCellTranslocationLymphToLung) and
                     e._moving_compartment == TBPulmonaryNetwork.T_CELL_ACTIVATED]
        self.assertEqual(len(tc_transl), 1)

        # T-cell death
        tn_death = [e for e in self.dynamics._events if isinstance(e, CellDeath) and
                    e._dying_compartment == TBPulmonaryNetwork.T_CELL_NAIVE]
        self.assertEqual(len(tn_death), 1)
        ta_death = [e for e in self.dynamics._events if isinstance(e, CellDeath) and
                    e._dying_compartment == TBPulmonaryNetwork.T_CELL_ACTIVATED]
        self.assertEqual(len(ta_death), 1)

    def test_configure_setUp_run(self):

        # Assign random param values
        params = {}

        params[TBPulmonaryNetwork.VENTILATION_SKEW] = 1.1
        params[TBPulmonaryNetwork.PERFUSION_SKEW] = 2.2
        params[TBPulmonaryNetwork.DRAINAGE_SKEW] = 3.3

        params[TBDynamics.IC_BAC_LOCATION] = 1
        params[TBDynamics.IC_BER_LOAD] = 1
        params[TBDynamics.IC_BED_LOAD] = 2

        # Import config file
        config = ConfigParser.ConfigParser()
        config.read('test_config.ini')

        # Assign parameters to the experiment
        for section in ['Event_Parameters', 'Network_Parameters', 'Initial_Conditions']:
            for k, v in config.items(section):
                # try to split the value by ",", if no "," exist, then it's just one value so use that value
                param = v.split(",")
                if len(param) == 1:
                    params[k] = float(v)
                else:
                    raise Exception("Invalid parameter format: {0}".format(v))

        self.dynamics.set_maximum_time(10)

        self.dynamics.set(params)
        self.dynamics.run()


if __name__ == '__main__':
    unittest.main()
