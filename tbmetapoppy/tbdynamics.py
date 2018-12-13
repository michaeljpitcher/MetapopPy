from metapoppy.dynamics import Dynamics
from .tbpulmonarynetwork import *
from tbmetapoppy.events import *


class TBDynamics(Dynamics):
    # Attribute seeding params
    VENTILATION_SKEW = 'ventilation_skew'
    PERFUSION_SKEW = 'perfusion_skew'
    DRAINAGE_SKEW = 'drainage_skew'

    # Event Parameter keys
    MACROPHAGE_CAPACITY = 'macrophage_capacity'
    SIGMOID_BIM_REPLICATION = 'sigmoid_intracellular_bacteria_replication'
    RP_REPLICATION_BER = 'replicating_extracellular_bacteria_replication_rate'
    RP_REPLICATION_BED = 'dormant_extracellular_bacteria_replication_rate'
    RP_REPLICATION_BIM = 'intracellular_bacteria_replication_rate'
    RP_CHANGE_TO_DORMANT = 'change_to_dormancy_rate'
    SIGMOID_CHANGE_TO_DORMANT = 'change_to_dormancy_sigmoid'
    HALF_SAT_CHANGE_TO_DORMANT = 'change_to_dormancy_halfsat'
    RP_CHANGE_TO_REPLICATING = 'change_to_replicating_rate'
    SIGMOID_CHANGE_TO_REPLICATING = 'change_to_replicating_sigmoid'
    HALF_SAT_CHANGE_TO_REPLICATING = 'change_to_replicating_halfsat'
    RP_BACTERIA_BLOOD_TRANSLOCATION = 'bacteria_blood_translocation_rate'

    RP_DCI_DEATH = 'immature_dendritic_cell_death_rate'
    RP_DCM_DEATH = 'mature_dendritic_cell_death_rate'
    RP_DCI_INGEST_BACTERIUM = 'dendritic_cell_ingest_bacterium_rate'
    HALF_SAT_DCI_INGEST_BACTERIUM = 'dendritic_cell_ingest_bacterium_halfsat'
    PROB_DC_INFECTION = 'dendritic_cell_infection_probability'
    RP_DCI_RECRUITMENT = 'dendritic_cell_recruitment_rate'
    RP_DCI_RECRUITMENT_ENHANCED = 'dendritic_cell_enhanced_recruitment_rate'
    HALF_SAT_DCI_RECRUITMENT_ENHANCED = 'dendritic_cell_enhanced_recruitment_half_sat'
    RP_DCM_TRANSLOCATION = 'dendritic_cell_translocation_rate'

    RP_MR_ACTIVATION_BY_TCELL = 'macrophage_activation_by_tcell_rate'
    HALF_SAT_MR_ACTIVATION_BY_TCELL = 'macrophage_activation_by_tcell_half_sat'
    RP_MR_ACTIVATION_BY_BACTERIA = 'macrophage_activation_by_bacteria_rate'
    HALF_SAT_MR_ACTIVATION_BY_BACTERIA = 'macrophage_activation_by_bacteria_half_sat'
    RP_MR_INGEST_BAC = 'resting_macrophage_ingest_bacterium_rate'
    HALF_SAT_MR_INGEST_BAC = 'resting_macrophage_ingest_bacterium_half_sat'
    PROB_MR_INFECTION = 'resting_macrophage_infection_probability'
    RP_MA_INGEST_BAC = 'activated_macrophage_ingest_bacterium_rate'
    HALF_SAT_MA_INGEST_BAC = 'activated_macrophage_ingest_bacterium_half_sat'
    PROB_MA_INFECTION = 'activated_macrophage_infection_probability'

    RP_MR_DEATH = 'resting_macrophage_death_rate'
    RP_MI_DEATH = 'infected_macrophage_death_rate'
    RP_MA_DEATH = 'activated_macrophage_death_rate'
    RP_MI_BURSTING = 'macrophage_bursting_rate'
    RP_TA_KILL_MI = 't_cell_destroys_macrophage_rate'
    HALF_SAT_TA_KILL_MI = 't_cell_destroys_macrophage_half_sat'
    RP_MR_RECRUIT_LUNG = 'macrophage_lung_recruitment_rate'
    RP_MR_RECRUIT_ENHANCE = 'macrophage_recruitment_enhanced_rate'
    HALF_SAT_MR_RECRUIT_ENHANCE = 'macrophage_recruitment_enhanced_half_sat'

    MI_TO_MA_CHEMOKINE_WEIGHT = 'macrophage_infected_to_activated_chemokine_weight'

    RP_MR_RECRUIT_LYMPH = 'macrophage_lymph_recruitment_rate'
    RP_MI_TRANSLOCATION = 'macrophage_translocation_rate'

    RP_TN_RECRUIT = 't_cell_recruitment_rate'
    RP_TCELL_ACTIVATION = 't_cell_activation_rate'
    HALF_SAT_TCELL_ACTIVATION = 't_cell_activation_half_sat'
    RP_TN_DEATH = 't_cell_naive_death_rate'
    RP_TA_DEATH = 't_cell_activated_death_rate'
    RP_TA_TRANSLOCATION = 't_cell_translocation_rate'

    IC_BER_LOAD = 'initial_bacterial_load_replicating'
    IC_BED_LOAD = 'initial_bacterial_load_dormant'

    def __init__(self, network_config):
        # Build network
        pulmonary_network = TBPulmonaryNetwork(network_config)
        Dynamics.__init__(self, pulmonary_network)

        self._perf_seed = {}

    def _create_events(self):
        """
        Create all TB related events, and assign relevant parameters to them
        :return:
        """
        events = []

        # Bacteria replication
        ber_replication = Replication(TBDynamics.RP_REPLICATION_BER, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING)
        bed_replication = Replication(TBDynamics.RP_REPLICATION_BED, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT)
        bim_replication = IntracellularBacterialReplication(TBDynamics.RP_REPLICATION_BIM,
                                                                  TBDynamics.SIGMOID_BIM_REPLICATION,
                                                                  TBDynamics.MACROPHAGE_CAPACITY)

        events += [ber_replication, bed_replication, bim_replication]

        # Bacterial change
        bed_to_ber = BacteriumChangeStateThroughOxygen(False, TBDynamics.RP_CHANGE_TO_REPLICATING,
                                                             TBDynamics.SIGMOID_CHANGE_TO_REPLICATING,
                                                             TBDynamics.HALF_SAT_CHANGE_TO_REPLICATING)
        ber_to_bed = BacteriumChangeStateThroughOxygen(True, TBDynamics.RP_CHANGE_TO_DORMANT,
                                                             TBDynamics.SIGMOID_CHANGE_TO_DORMANT,
                                                             TBDynamics.HALF_SAT_CHANGE_TO_DORMANT)
        events += [bed_to_ber, ber_to_bed]

        # bacterial translocation
        bed_translocation = CellTranslocationToLung(TBDynamics.RP_BACTERIA_BLOOD_TRANSLOCATION,
                                                    TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT)
        events.append(bed_translocation)

        # Dendritic cell recruitment
        dc_recruit_lung = StandardCellRecruitmentLung(TBDynamics.RP_DCI_RECRUITMENT, TBPulmonaryNetwork.DENDRITIC_CELL_IMMATURE)
        dc_recruit_lung_enh = EnhancedCellRecruitmentLung(TBDynamics.RP_DCI_RECRUITMENT_ENHANCED,
                                                          TBPulmonaryNetwork.DENDRITIC_CELL_IMMATURE,
                                                                     TBDynamics.HALF_SAT_DCI_RECRUITMENT_ENHANCED,
                                                                     TBDynamics.MI_TO_MA_CHEMOKINE_WEIGHT)
        events += [dc_recruit_lung, dc_recruit_lung_enh]

        # Dendritic cell death
        dci_death = CellDeath(TBDynamics.RP_DCI_DEATH, TBPulmonaryNetwork.DENDRITIC_CELL_IMMATURE)
        dcm_death = CellDeath(TBDynamics.RP_DCM_DEATH, TBPulmonaryNetwork.DENDRITIC_CELL_MATURE)
        events += [dci_death, dcm_death]

        # Dendritic Cell maturation
        dc_maturation = CellIngestBacterium(TBDynamics.RP_DCI_INGEST_BACTERIUM, TBPulmonaryNetwork.DENDRITIC_CELL_IMMATURE,
                                            TBDynamics.HALF_SAT_DCI_INGEST_BACTERIUM, TBDynamics.PROB_DC_INFECTION)
        events.append(dc_maturation)

        # Dendritic cell translocation
        dc_translocation = CellTranslocationToLymph(TBDynamics.RP_DCM_TRANSLOCATION, TBPulmonaryNetwork.DENDRITIC_CELL_MATURE)
        events.append(dc_translocation)

        # Macrophage recruitment
        mr_recruit_lung = StandardCellRecruitmentLung(TBDynamics.RP_MR_RECRUIT_LUNG, TBPulmonaryNetwork.MACROPHAGE_RESTING)
        mr_recruit_lymph = StandardCellRecruitmentLymph(TBDynamics.RP_MR_RECRUIT_LYMPH, TBPulmonaryNetwork.MACROPHAGE_RESTING)
        mr_recruit_lung_enhanced = EnhancedCellRecruitmentLung(TBDynamics.RP_MR_RECRUIT_ENHANCE,
                                                               TBPulmonaryNetwork.MACROPHAGE_RESTING,
                                                               TBDynamics.HALF_SAT_MR_RECRUIT_ENHANCE,
                                                               TBDynamics.MI_TO_MA_CHEMOKINE_WEIGHT)
        mr_recruit_lymph_enhanced = EnhancedCellRecruitmentLymph(TBDynamics.RP_MR_RECRUIT_ENHANCE,
                                                                 TBPulmonaryNetwork.MACROPHAGE_RESTING,
                                                               TBDynamics.HALF_SAT_MR_RECRUIT_ENHANCE,
                                                               TBDynamics.MI_TO_MA_CHEMOKINE_WEIGHT)
        
        events += [mr_recruit_lung, mr_recruit_lymph, mr_recruit_lung_enhanced, mr_recruit_lymph_enhanced]

        # Macrophage activation
        mr_activation_bac = CellActivation(TBDynamics.RP_MR_ACTIVATION_BY_BACTERIA,
                                           TBDynamics.HALF_SAT_MR_ACTIVATION_BY_BACTERIA, TBPulmonaryNetwork.MACROPHAGE_RESTING,
                                           TBPulmonaryNetwork.EXTRACELLULAR_BACTERIA)
        mr_activation_ta = CellActivation(TBDynamics.RP_MR_ACTIVATION_BY_TCELL,
                                          TBDynamics.HALF_SAT_MR_ACTIVATION_BY_TCELL,
                                          TBPulmonaryNetwork.MACROPHAGE_RESTING, [TBPulmonaryNetwork.T_CELL_ACTIVATED])
        events += [mr_activation_bac, mr_activation_ta]

        # Macrophage death
        mr_death = CellDeath(TBDynamics.RP_MR_DEATH, TBPulmonaryNetwork.MACROPHAGE_RESTING)
        ma_death = CellDeath(TBDynamics.RP_MA_DEATH, TBPulmonaryNetwork.MACROPHAGE_ACTIVATED)
        mi_death = CellDeath(TBDynamics.RP_MI_DEATH, TBPulmonaryNetwork.MACROPHAGE_INFECTED)
        mi_burst = MacrophageBursting(TBDynamics.RP_MI_BURSTING, TBDynamics.SIGMOID_BIM_REPLICATION,
                                      TBDynamics.MACROPHAGE_CAPACITY)

        ta_kills_mi = TCellDestroysMacrophage(TBDynamics.RP_TA_KILL_MI, TBDynamics.HALF_SAT_TA_KILL_MI)
        events += [mr_death, ma_death, mi_death, mi_burst, ta_kills_mi]

        # Macrophage ingest bacterium
        mr_ingest_bac = CellIngestBacterium(TBDynamics.RP_MR_INGEST_BAC, TBPulmonaryNetwork.MACROPHAGE_RESTING,
                                            TBDynamics.HALF_SAT_MR_INGEST_BAC,
                                            TBDynamics.PROB_MR_INFECTION)
        ma_ingest_bac = CellIngestBacterium(TBDynamics.RP_MA_INGEST_BAC, TBPulmonaryNetwork.MACROPHAGE_ACTIVATED,
                                            TBDynamics.HALF_SAT_MA_INGEST_BAC,
                                            TBDynamics.PROB_MA_INFECTION)
        events += [mr_ingest_bac, ma_ingest_bac]

        # Macrophage translocation
        mi_translocation = CellTranslocationToLymph(TBDynamics.RP_MI_TRANSLOCATION, TBPulmonaryNetwork.MACROPHAGE_INFECTED)
        events.append(mi_translocation)

        # T-cell recruitment
        tn_recruit = StandardCellRecruitmentLymph(TBDynamics.RP_TN_RECRUIT, TBPulmonaryNetwork.T_CELL_NAIVE)
        events.append(tn_recruit)

        # T-cell activation
        tn_activation = CellActivation(TBDynamics.RP_TCELL_ACTIVATION, TBDynamics.HALF_SAT_TCELL_ACTIVATION,
                                       TBPulmonaryNetwork.T_CELL_NAIVE, [TBPulmonaryNetwork.DENDRITIC_CELL_MATURE, TBPulmonaryNetwork.MACROPHAGE_INFECTED])
        events.append(tn_activation)

        # T-cell translocation
        ta_translocation = TCellTranslocationToLungByInfection(TBDynamics.RP_TA_TRANSLOCATION)
        events.append(ta_translocation)

        # T-cell death
        tn_death = CellDeath(TBDynamics.RP_TN_DEATH, TBPulmonaryNetwork.T_CELL_NAIVE)
        ta_death = CellDeath(TBDynamics.RP_TA_DEATH, TBPulmonaryNetwork.T_CELL_ACTIVATED)
        events += [tn_death, ta_death]

        return events

    def _get_patch_seeding(self, params):
        # Get patch attributes
        patch_seeding = self._network.get_pulmonary_att_seeding(params[TBDynamics.VENTILATION_SKEW],
                                                                    params[TBDynamics.PERFUSION_SKEW],
                                                                    params[TBDynamics.DRAINAGE_SKEW])

        # Compartment seeding
        mac_recruit_lung = params[TBDynamics.RP_MR_RECRUIT_LUNG]
        mac_recruit_lymph = params[TBDynamics.RP_MR_RECRUIT_LYMPH]
        mac_death = params[TBDynamics.RP_MR_DEATH]
        dc_recruit = params[TBDynamics.RP_DCI_RECRUITMENT]
        dc_death = params[TBDynamics.RP_DCI_DEATH]
        tn_recruit = params[TBDynamics.RP_TN_RECRUIT]
        tn_death = params[TBDynamics.RP_TN_DEATH]
        for n, v in patch_seeding.iteritems():
            self._perf_seed[n] = v[TypedNetwork.ATTRIBUTES][TBPulmonaryNetwork.PERFUSION]
            v[TypedNetwork.COMPARTMENTS] = {TBPulmonaryNetwork.MACROPHAGE_RESTING: int(round(self._perf_seed[n] *
                                                                          (mac_recruit_lung / mac_death))),
                                            TBPulmonaryNetwork.DENDRITIC_CELL_IMMATURE: int(round(self._perf_seed[n] *
                                                                               (dc_recruit / dc_death)))}
        lymph_seed = {TBPulmonaryNetwork.MACROPHAGE_RESTING: int(round(mac_recruit_lymph / mac_death)),
                      TBPulmonaryNetwork.T_CELL_NAIVE: int(round(tn_recruit / tn_death))}
        patch_seeding[TBPulmonaryNetwork.LYMPH_PATCH] = {TypedNetwork.COMPARTMENTS: lymph_seed}

        # Bacteria
        # Ventilation based - assumes sum of ventilation values = 1.0
        # r = numpy.random.random()
        # count = 0
        # for p in self._network.get_patches_by_type(PulmonaryNetwork.ALVEOLAR_PATCH):
        #     count += patch_seeding[p][TypedNetwork.ATTRIBUTES][PulmonaryNetwork.VENTILATION]
        #     if count >= r:
        #         patch_seeding[p][TypedNetwork.COMPARTMENTS][BACTERIUM_EXTRACELLULAR_REPLICATING] = \
        #             params[TBDynamics.IC_BER_LOAD]
        #         patch_seeding[p][TypedNetwork.COMPARTMENTS][BACTERIUM_EXTRACELLULAR_DORMANT] = \
        #             params[TBDynamics.IC_BED_LOAD]
        #         break
        # TODO - debug code
        patch_seeding[9][TypedNetwork.COMPARTMENTS][TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING] = \
            params[TBDynamics.IC_BER_LOAD]
        patch_seeding[9][TypedNetwork.COMPARTMENTS][TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT] = \
            params[TBDynamics.IC_BED_LOAD]
        return patch_seeding

    def _build_network(self, params):
        pass

    def _get_edge_seeding(self, params):
        seeding = {}
        for n, v in self._network.get_patches_by_type(TBPulmonaryNetwork.ALVEOLAR_PATCH, data=True):
            seeding[(n, TBPulmonaryNetwork.LYMPH_PATCH)] = {TBPulmonaryNetwork.PERFUSION: self._perf_seed[n]}
        return seeding

    def _patch_is_active(self, patch_id):
        """
        Determine if the given patch is active (from the network). Alveolar patches only become active when they contain
        bacteria.
        :param patch_id:
        :return:
        """
        return patch_id == TBPulmonaryNetwork.LYMPH_PATCH or \
               sum([self._network.get_compartment_value(patch_id, n) for n in TBPulmonaryNetwork.BACTERIA]) > 0
