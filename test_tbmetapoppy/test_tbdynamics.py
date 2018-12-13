import unittest
from tbmetapoppy import *


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
        bed_transl = [e for e in self.dynamics._events if isinstance(e, CellTranslocationToLung) and
                      e._cell_type == TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT]
        self.assertEqual(len(bed_transl), 1)

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
        dc_transl = [e for e in self.dynamics._events if isinstance(e, CellTranslocationToLymph) and
                     e._cell_type == TBPulmonaryNetwork.DENDRITIC_CELL_MATURE]
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
        mi_nat_death = [e for e in self.dynamics._events if e.__class__ == CellDeath and
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
        mac_transl = [e for e in self.dynamics._events if isinstance(e, CellTranslocationToLymph) and
                      e._cell_type == TBPulmonaryNetwork.MACROPHAGE_INFECTED]
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
        tc_transl = [e for e in self.dynamics._events if isinstance(e, CellTranslocationToLung) and
                     e._cell_type == TBPulmonaryNetwork.T_CELL_ACTIVATED]
        self.assertEqual(len(tc_transl), 1)

        # T-cell death
        tn_death = [e for e in self.dynamics._events if isinstance(e, CellDeath) and
                    e._dying_compartment == TBPulmonaryNetwork.T_CELL_NAIVE]
        self.assertEqual(len(tn_death), 1)
        ta_death = [e for e in self.dynamics._events if isinstance(e, CellDeath) and
                    e._dying_compartment == TBPulmonaryNetwork.T_CELL_ACTIVATED]
        self.assertEqual(len(ta_death), 1)

    def test_configure_setUp_run(self):
        params = {}
        params[TBDynamics.RP_REPLICATION_BER] = 0.01
        params[TBDynamics.RP_REPLICATION_BED] = 0.02
        params[TBDynamics.RP_REPLICATION_BIM] = 0.03
        params[TBDynamics.SIGMOID_BIM_REPLICATION] = 0.04
        params[TBDynamics.MACROPHAGE_CAPACITY] = 0.05
        params[TBDynamics.RP_CHANGE_TO_REPLICATING] = 0.06
        params[TBDynamics.SIGMOID_CHANGE_TO_REPLICATING] = 2
        params[TBDynamics.HALF_SAT_CHANGE_TO_REPLICATING] = 0.07
        params[TBDynamics.RP_CHANGE_TO_DORMANT] = 0.08
        params[TBDynamics.SIGMOID_CHANGE_TO_DORMANT] = -2
        params[TBDynamics.HALF_SAT_CHANGE_TO_DORMANT] = 0.09
        params[TBDynamics.RP_MR_INGEST_BAC] = 0.10
        params[TBDynamics.HALF_SAT_MR_INGEST_BAC] = 111
        params[TBDynamics.PROB_MR_INFECTION] = 0.5
        params[TBDynamics.RP_MA_INGEST_BAC] = 0.11
        params[TBDynamics.HALF_SAT_MA_INGEST_BAC] = 112
        params[TBDynamics.PROB_MA_INFECTION] = 0.0
        params[TBDynamics.RP_BACTERIA_BLOOD_TRANSLOCATION] = 0.12
        params[TBDynamics.RP_DCI_RECRUITMENT] = 130.0
        params[TBDynamics.RP_DCI_DEATH] = 0.14
        params[TBDynamics.RP_DCM_DEATH] = 0.15
        params[TBDynamics.RP_DCI_INGEST_BACTERIUM] = 0.16
        params[TBDynamics.HALF_SAT_DCI_INGEST_BACTERIUM] = 0.17
        params[TBDynamics.PROB_DC_INFECTION] = 1.0
        params[TBDynamics.RP_DCM_TRANSLOCATION] = 0.18
        params[TBDynamics.RP_MR_RECRUIT_LUNG] = 190.0
        params[TBDynamics.RP_MR_RECRUIT_LYMPH] = 200.0
        params[TBDynamics.RP_MR_ACTIVATION_BY_TCELL] = 0.21
        params[TBDynamics.HALF_SAT_MR_ACTIVATION_BY_TCELL] = 0.22
        params[TBDynamics.RP_MR_ACTIVATION_BY_BACTERIA] = 0.23
        params[TBDynamics.HALF_SAT_MR_ACTIVATION_BY_BACTERIA] = 0.24
        params[TBDynamics.RP_MR_DEATH] = 0.25
        params[TBDynamics.RP_MA_DEATH] = 0.26
        params[TBDynamics.RP_MI_DEATH] = 0.27
        params[TBDynamics.RP_MI_BURSTING] = 0.28
        params[TBDynamics.RP_TA_KILL_MI] = 0.29
        params[TBDynamics.HALF_SAT_TA_KILL_MI] = 0.30
        params[TBDynamics.RP_MI_TRANSLOCATION] = 0.33
        params[TBDynamics.RP_TN_RECRUIT] = 340.0
        params[TBDynamics.RP_TCELL_ACTIVATION] = 0.35
        params[TBDynamics.HALF_SAT_TCELL_ACTIVATION] = 0.36
        params[TBDynamics.RP_TA_TRANSLOCATION] = 0.37
        params[TBDynamics.RP_TN_DEATH] = 0.38
        params[TBDynamics.RP_TA_DEATH] = 0.39

        params[TBDynamics.RP_MR_RECRUIT_ENHANCE] = 0.42
        params[TBDynamics.HALF_SAT_MR_RECRUIT_ENHANCE] = 0.43
        params[TBDynamics.MI_TO_MA_CHEMOKINE_WEIGHT] = 0.5

        params[TBDynamics.RP_DCI_RECRUITMENT_ENHANCED] = 0.44
        params[TBDynamics.HALF_SAT_DCI_RECRUITMENT_ENHANCED] = 0.43

        params[TBDynamics.VENTILATION_SKEW] = 1.1
        params[TBDynamics.PERFUSION_SKEW] = 2.2
        params[TBDynamics.DRAINAGE_SKEW] = 3.3

        params[TBDynamics.IC_BER_LOAD] = 1
        params[TBDynamics.IC_BED_LOAD] = 2

        self.dynamics.configure(params)
        # TODO = Check values stored correctly
        # Bacterial replication

        ber_rep = next(e for e in self.dynamics._events if isinstance(e, Replication) and
                   e._cell_type == TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING)
        self.assertEqual(ber_rep._reaction_parameter, params[TBDynamics.RP_REPLICATION_BER])

        bed_rep = next(e for e in self.dynamics._events if isinstance(e, Replication) and
                   e._cell_type == TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT)
        self.assertEqual(bed_rep._reaction_parameter, params[TBDynamics.RP_REPLICATION_BED])

        bim_rep = next(e for e in self.dynamics._events if isinstance(e, IntracellularBacterialReplication))
        self.assertEqual(bim_rep._reaction_parameter, params[TBDynamics.RP_REPLICATION_BIM])
        self.assertEqual(bim_rep._parameters[TBDynamics.SIGMOID_BIM_REPLICATION],
                         params[TBDynamics.SIGMOID_BIM_REPLICATION])
        self.assertEqual(bim_rep._parameters[TBDynamics.MACROPHAGE_CAPACITY],
                         params[TBDynamics.MACROPHAGE_CAPACITY])

        # Bacterial state change
        bed_to_ber = next(e for e in self.dynamics._events if isinstance(e, BacteriumChangeStateThroughOxygen)
                          and e._compartment_from == TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT)
        self.assertEqual(bed_to_ber._reaction_parameter, params[TBDynamics.RP_CHANGE_TO_REPLICATING])
        self.assertEqual(bed_to_ber._parameters[TBDynamics.SIGMOID_CHANGE_TO_REPLICATING],
                         params[TBDynamics.SIGMOID_CHANGE_TO_REPLICATING])
        self.assertEqual(bed_to_ber._parameters[TBDynamics.HALF_SAT_CHANGE_TO_REPLICATING],
                         params[TBDynamics.HALF_SAT_CHANGE_TO_REPLICATING])
        ber_to_bed = next(e for e in self.dynamics._events if isinstance(e, BacteriumChangeStateThroughOxygen)
                          and e._compartment_from == TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING)
        self.assertEqual(ber_to_bed._reaction_parameter, params[TBDynamics.RP_CHANGE_TO_DORMANT])
        self.assertEqual(ber_to_bed._parameters[TBDynamics.SIGMOID_CHANGE_TO_DORMANT],
                         params[TBDynamics.SIGMOID_CHANGE_TO_DORMANT])
        self.assertEqual(ber_to_bed._parameters[TBDynamics.HALF_SAT_CHANGE_TO_DORMANT],
                         params[TBDynamics.HALF_SAT_CHANGE_TO_DORMANT])

        # Bacterial translocation
        bed_transl = next(e for e in self.dynamics._events if isinstance(e, CellTranslocationToLung) and
                          e._cell_type == TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT)
        self.assertEqual(bed_transl._reaction_parameter, params[TBDynamics.RP_BACTERIA_BLOOD_TRANSLOCATION])

        # DC recruitment
        dc_rec_lung = next(e for e in self.dynamics._events if isinstance(e, StandardCellRecruitmentLung) and
                           e._cell_type == TBPulmonaryNetwork.DENDRITIC_CELL_IMMATURE)
        self.assertEqual(dc_rec_lung._reaction_parameter, params[TBDynamics.RP_DCI_RECRUITMENT])

        # DC infection
        dc_infect = next(e for e in self.dynamics._events if isinstance(e, CellIngestBacterium) and
                     e._cell_type == TBPulmonaryNetwork.DENDRITIC_CELL_IMMATURE)
        self.assertEqual(dc_infect._reaction_parameter, params[TBDynamics.RP_DCI_INGEST_BACTERIUM])
        self.assertEqual(dc_infect._parameters[TBDynamics.HALF_SAT_DCI_INGEST_BACTERIUM],
                         params[TBDynamics.HALF_SAT_DCI_INGEST_BACTERIUM])
        self.assertEqual(dc_infect._parameters[TBDynamics.PROB_DC_INFECTION],
                         params[TBDynamics.PROB_DC_INFECTION])

        # DC translocation
        dc_transl = next(e for e in self.dynamics._events if isinstance(e, CellTranslocationToLymph) and
                     e._cell_type == TBPulmonaryNetwork.DENDRITIC_CELL_MATURE)
        self.assertEqual(dc_transl._reaction_parameter, params[TBDynamics.RP_DCM_TRANSLOCATION])

        # DC death
        dci_death = next(e for e in self.dynamics._events if isinstance(e, CellDeath) and
                         e._dying_compartment == TBPulmonaryNetwork.DENDRITIC_CELL_IMMATURE)
        self.assertEqual(dci_death._reaction_parameter, params[TBDynamics.RP_DCI_DEATH])
        dcm_death = next(e for e in self.dynamics._events if isinstance(e, CellDeath) and
                         e._dying_compartment == TBPulmonaryNetwork.DENDRITIC_CELL_MATURE)
        self.assertEqual(dcm_death._reaction_parameter, params[TBDynamics.RP_DCM_DEATH])

        # Macrophage recruitment
        mac_rec_lung = next(e for e in self.dynamics._events if isinstance(e, StandardCellRecruitmentLung) and
                            e._cell_type == TBPulmonaryNetwork.MACROPHAGE_RESTING)
        self.assertEqual(mac_rec_lung._reaction_parameter, params[TBDynamics.RP_MR_RECRUIT_LUNG])
        mac_rec_lymph = next(e for e in self.dynamics._events if isinstance(e, StandardCellRecruitmentLymph) and
                             e._cell_type == TBPulmonaryNetwork.MACROPHAGE_RESTING)
        self.assertEqual(mac_rec_lymph._reaction_parameter, params[TBDynamics.RP_MR_RECRUIT_LYMPH])

        mac_rec_lung_enh = next(e for e in self.dynamics._events if isinstance(e, EnhancedCellRecruitmentLung) and
                                   e._cell_type == TBPulmonaryNetwork.MACROPHAGE_RESTING)
        self.assertEqual(mac_rec_lung_enh._reaction_parameter, params[TBDynamics.RP_MR_RECRUIT_ENHANCE])
        self.assertEqual(mac_rec_lung_enh._parameters[TBDynamics.HALF_SAT_MR_RECRUIT_ENHANCE],
                         params[TBDynamics.HALF_SAT_MR_RECRUIT_ENHANCE])
        self.assertEqual(mac_rec_lung_enh._parameters[TBDynamics.MI_TO_MA_CHEMOKINE_WEIGHT],
                         params[TBDynamics.MI_TO_MA_CHEMOKINE_WEIGHT])


        mac_rec_lym_enh = next(e for e in self.dynamics._events if isinstance(e, EnhancedCellRecruitmentLymph) and
                                  e._cell_type == TBPulmonaryNetwork.MACROPHAGE_RESTING)
        self.assertEqual(mac_rec_lym_enh._reaction_parameter, params[TBDynamics.RP_MR_RECRUIT_ENHANCE])
        self.assertEqual(mac_rec_lym_enh._parameters[TBDynamics.HALF_SAT_MR_RECRUIT_ENHANCE],
                         params[TBDynamics.HALF_SAT_MR_RECRUIT_ENHANCE])
        self.assertEqual(mac_rec_lym_enh._parameters[TBDynamics.MI_TO_MA_CHEMOKINE_WEIGHT],
                         params[TBDynamics.MI_TO_MA_CHEMOKINE_WEIGHT])

        # Macrophage activation
        mac_act_tc = next(e for e in self.dynamics._events if isinstance(e, CellActivation) and
                      e._resting_cell == TBPulmonaryNetwork.MACROPHAGE_RESTING and e._triggers == [TBPulmonaryNetwork.T_CELL_ACTIVATED])
        self.assertEqual(mac_act_tc._reaction_parameter, params[TBDynamics.RP_MR_ACTIVATION_BY_TCELL])
        self.assertEqual(mac_act_tc._parameters[TBDynamics.HALF_SAT_MR_ACTIVATION_BY_TCELL],
                         params[TBDynamics.HALF_SAT_MR_ACTIVATION_BY_TCELL])
        mac_act_bac = next(e for e in self.dynamics._events if isinstance(e, CellActivation) and
                       e._resting_cell == TBPulmonaryNetwork.MACROPHAGE_RESTING and
                       (e._triggers == [TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING]
                        or e._triggers == [TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT]))
        self.assertEqual(mac_act_bac._reaction_parameter, params[TBDynamics.RP_MR_ACTIVATION_BY_BACTERIA])
        self.assertEqual(mac_act_bac._parameters[TBDynamics.HALF_SAT_MR_ACTIVATION_BY_BACTERIA],
                         params[TBDynamics.HALF_SAT_MR_ACTIVATION_BY_BACTERIA])

        # Macrophage death
        mr_death = next(e for e in self.dynamics._events if isinstance(e, CellDeath)
                        and e._dying_compartment == TBPulmonaryNetwork.MACROPHAGE_RESTING)
        self.assertEqual(mr_death._reaction_parameter, params[TBDynamics.RP_MR_DEATH])
        ma_death = next(e for e in self.dynamics._events if isinstance(e, CellDeath)
                        and e._dying_compartment == TBPulmonaryNetwork.MACROPHAGE_ACTIVATED)
        self.assertEqual(ma_death._reaction_parameter, params[TBDynamics.RP_MA_DEATH])

        # Not isinstance since bursting and t-cell death are instances (sub-classes)
        mi_nat_death = next(e for e in self.dynamics._events if e.__class__ == CellDeath and
                        e._dying_compartment == TBPulmonaryNetwork.MACROPHAGE_INFECTED)
        self.assertEqual(mi_nat_death._reaction_parameter, params[TBDynamics.RP_MI_DEATH])
        mi_burst = next(e for e in self.dynamics._events if isinstance(e, MacrophageBursting))
        self.assertEqual(mi_burst._reaction_parameter, params[TBDynamics.RP_MI_BURSTING])
        self.assertEqual(mi_burst._parameters[TBDynamics.SIGMOID_BIM_REPLICATION],
                         params[TBDynamics.SIGMOID_BIM_REPLICATION])
        self.assertEqual(mi_burst._parameters[TBDynamics.MACROPHAGE_CAPACITY],
                         params[TBDynamics.MACROPHAGE_CAPACITY])
        mi_death_tcell = next(e for e in self.dynamics._events if isinstance(e, TCellDestroysMacrophage))
        self.assertEqual(mi_death_tcell._reaction_parameter, params[TBDynamics.RP_TA_KILL_MI])
        self.assertEqual(mi_death_tcell._parameters[TBDynamics.HALF_SAT_TA_KILL_MI],
                         params[TBDynamics.HALF_SAT_TA_KILL_MI])

        # Macrophage ingest
        mr_ingest_bac = next(e for e in self.dynamics._events if isinstance(e, CellIngestBacterium) and
                             e._cell_type == TBPulmonaryNetwork.MACROPHAGE_RESTING)
        self.assertEqual(mr_ingest_bac._reaction_parameter, params[TBDynamics.RP_MR_INGEST_BAC])
        self.assertEqual(mr_ingest_bac._parameters[TBDynamics.HALF_SAT_MR_INGEST_BAC],
                         params[TBDynamics.HALF_SAT_MR_INGEST_BAC])
        self.assertEqual(mr_ingest_bac._parameters[TBDynamics.PROB_MR_INFECTION],
                         params[TBDynamics.PROB_MR_INFECTION])
        ma_ingest_bac = next(e for e in self.dynamics._events if isinstance(e, CellIngestBacterium) and
                             e._cell_type == TBPulmonaryNetwork.MACROPHAGE_ACTIVATED)
        self.assertEqual(ma_ingest_bac._reaction_parameter, params[TBDynamics.RP_MA_INGEST_BAC])
        self.assertEqual(ma_ingest_bac._parameters[TBDynamics.HALF_SAT_MA_INGEST_BAC],
                         params[TBDynamics.HALF_SAT_MA_INGEST_BAC])
        self.assertEqual(ma_ingest_bac._parameters[TBDynamics.PROB_MA_INFECTION],
                         params[TBDynamics.PROB_MA_INFECTION])

        # Macrophage translocation
        mac_transl = next(e for e in self.dynamics._events if isinstance(e, CellTranslocationToLymph) and
                      e._cell_type == TBPulmonaryNetwork.MACROPHAGE_INFECTED)
        self.assertEqual(mac_transl._reaction_parameter, params[TBDynamics.RP_MI_TRANSLOCATION])

        #  T-cell recruitment
        tc_recruit = next(e for e in self.dynamics._events if isinstance(e, StandardCellRecruitmentLymph) and
                          e._cell_type == TBPulmonaryNetwork.T_CELL_NAIVE)
        self.assertEqual(tc_recruit._reaction_parameter, params[TBDynamics.RP_TN_RECRUIT])

        # T-cell activation
        tc_act = next(e for e in self.dynamics._events if isinstance(e, CellActivation) and
                  e._resting_cell == TBPulmonaryNetwork.T_CELL_NAIVE)
        self.assertEqual(tc_act._reaction_parameter, params[TBDynamics.RP_TCELL_ACTIVATION])
        self.assertEqual(tc_act._parameters[TBDynamics.HALF_SAT_TCELL_ACTIVATION],
                         params[TBDynamics.HALF_SAT_TCELL_ACTIVATION])

        # T-cell translocation
        tc_transl = next(e for e in self.dynamics._events if isinstance(e, CellTranslocationToLung) and
                     e._cell_type == TBPulmonaryNetwork.T_CELL_ACTIVATED)
        self.assertEqual(tc_transl._reaction_parameter, params[TBDynamics.RP_TA_TRANSLOCATION])

        # T-cell death
        tn_death = next(e for e in self.dynamics._events if isinstance(e, CellDeath) and
                        e._dying_compartment == TBPulmonaryNetwork.T_CELL_NAIVE)
        self.assertEqual(tn_death._reaction_parameter, params[TBDynamics.RP_TN_DEATH])
        ta_death = next(e for e in self.dynamics._events if isinstance(e, CellDeath) and
                        e._dying_compartment == TBPulmonaryNetwork.T_CELL_ACTIVATED)
        self.assertEqual(ta_death._reaction_parameter, params[TBDynamics.RP_TA_DEATH])

        # Set Up
        self.dynamics.setUp(params)

        # Test seed by rates
        for p, v in self.dynamics.network().nodes(data=True):
            if v[TypedNetwork.PATCH_TYPE] == TBPulmonaryNetwork.ALVEOLAR_PATCH:
                # TODO test pul atts
                pass
                # Test comp seed
                self.assertEqual(v[TBPulmonaryNetwork.COMPARTMENTS][TBPulmonaryNetwork.MACROPHAGE_RESTING],
                                 int(round(v[TBPulmonaryNetwork.ATTRIBUTES][TBPulmonaryNetwork.PERFUSION] *
                                           (params[TBDynamics.RP_MR_RECRUIT_LUNG]) / params[TBDynamics.RP_MR_DEATH])))
                self.assertEqual(v[TBPulmonaryNetwork.COMPARTMENTS][TBPulmonaryNetwork.DENDRITIC_CELL_IMMATURE],
                                 int(round(v[TBPulmonaryNetwork.ATTRIBUTES][TBPulmonaryNetwork.PERFUSION] *
                                           (params[TBDynamics.RP_DCI_RECRUITMENT]) / params[TBDynamics.RP_DCI_DEATH])))
            elif v[TypedNetwork.PATCH_TYPE] == TBPulmonaryNetwork.LYMPH_PATCH:
                self.assertEqual(v[TBPulmonaryNetwork.COMPARTMENTS][TBPulmonaryNetwork.T_CELL_NAIVE],
                                 int(round(params[TBDynamics.RP_TN_RECRUIT] / params[TBDynamics.RP_TN_DEATH])))
            else:
                raise Exception

        # Do
        self.dynamics.set_maximum_time(10)
        self.dynamics.do(params)


if __name__ == '__main__':
    unittest.main()
