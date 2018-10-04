from metapoppy.dynamics import Dynamics
from .pulmonarynetwork import *
from tbmetapoppy.events import *
from tbcompartments import *


class TBDynamics(Dynamics):
    # TODO - enhanced recruitment rates

    # Event Parameter keys
    MACROPHAGE_CAPACITY = 'macrophage_capacity'
    SIGMOID_BIM_REPLICATION = 'sigmoid_intracellular_bacteria_replication'
    RP_REPLICATION_BER = 'replicating_extracellular_bacteria_replication_rate'
    RP_REPLICATION_BED = 'dormant_extracellular_bacteria_replication_rate'
    RP_REPLICATION_BIM = 'intracellular_bacteria_replication_rate'
    RP_MR_DESTROY_BACTERIA = 'regular_macrophage_bacterial_destruction_rate'
    HALF_SAT_MR_DESTROY_BACTERIA = 'regular_macrophage_bacterial_destruction_half_sat'
    RP_MA_DESTROY_BACTERIA = 'activated_macrophage_bacterial_destruction_rate'
    HALF_SAT_MA_DESTROY_BACTERIA = 'activated_macrophage_bacterial_destruction_half_sat'
    RP_CHANGE_TO_DORMANT = 'change_to_dormancy_rate'
    SIGMOID_CHANGE_TO_DORMANT = 'change_to_dormancy_sigmoid'
    HALF_SAT_CHANGE_TO_DORMANT = 'change_to_dormancy_halfsat'
    RP_CHANGE_TO_REPLICATING = 'change_to_replicating_rate'
    SIGMOID_CHANGE_TO_REPLICATING = 'change_to_replicating_sigmoid'
    HALF_SAT_CHANGE_TO_REPLICATING = 'change_to_replicating_halfsat'
    RP_BACTERIA_BLOOD_TRANSLOCATION = 'bacteria_blood_translocation_rate'

    RP_DCI_DEATH = 'immature_dendritic_cell_death_rate'
    RP_DCM_DEATH = 'mature_dendritic_cell_death_rate'
    RP_DCI_MATURATION = 'dendritic_cell_maturation_rate'
    HALF_SAT_DCI_MATURATION = 'dendritic_cell_maturation_halfsat'
    RP_DCI_RECRUITMENT = 'dendritic_cell_recruitment_rate'
    RP_DCM_TRANSLOCATION = 'dendrtic_cell_translocation_rate'

    RP_MR_ACTIVATION_BY_TCELL = 'macrophage_activation_by_tcell_rate'
    HALF_SAT_MR_ACTIVATION_BY_TCELL = 'macrophage_activation_by_tcell_half_sat'
    RP_MR_ACTIVATION_BY_BACTERIA = 'macrophage_activation_by_bacteria_rate'
    HALF_SAT_MR_ACTIVATION_BY_BACTERIA = 'macrophage_activation_by_bacteria_half_sat'
    RP_MR_INFECTION = 'macrophage_infection_rate'
    HALF_SAT_MR_INFECTION = 'macrophage_infection_half_sat'
    RP_MR_DEATH = 'resting_macrophage_death_rate'
    RP_MI_DEATH = 'infected_macrophage_death_rate'
    RP_MA_DEATH = 'activated_macrophage_death_rate'
    RP_MI_BURSTING = 'macrophage_bursting_rate'
    RP_TA_KILL_MI = 't_cell_destroys_macrophage_rate'
    HALF_SAT_TA_KILL_MI = 't_cell_destroys_macrophage_half_sat'
    RP_MR_RECRUIT_LUNG = 'macrophage_lung_recruitment_rate'
    RP_MR_RECRUIT_LYMPH = 'macrophage_lymph_recruitment_rate'
    RP_MI_TRANSLOCATION = 'macrophage_translocation_rate'

    RP_TN_RECRUIT = 't_cell_recruitment_rate'
    RP_TCELL_ACTIVATION = 't_cell_activation_rate'
    HALF_SAT_TCELL_ACTIVATION = 't_cell_activation_half_sat'
    RP_TN_DEATH = 't_cell_naive_death_rate'
    RP_TA_DEATH = 't_cell_activated_death_rate'
    RP_TA_TRANSLOCATION = 't_cell_translocation_rate'

    IC_BACTERIA_LOADS = 'initial_bacterial_loads'

    def __init__(self, network_config):
        # Build network
        pulmonary_network = PulmonaryNetwork(network_config, TB_COMPARTMENTS)
        Dynamics.__init__(self, pulmonary_network)

    def _create_events(self):
        events = []

        # Bacteria replication
        ber_replication = Replication(TBDynamics.RP_REPLICATION_BER, BACTERIUM_EXTRACELLULAR_REPLICATING)
        bed_replication = Replication(TBDynamics.RP_REPLICATION_BED, BACTERIUM_EXTRACELLULAR_DORMANT)
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

        # Bacterial destruction
        mr_kills_bac = MacrophageDestroysBacterium(TBDynamics.RP_MR_DESTROY_BACTERIA, MACROPHAGE_RESTING,
                                                   TBDynamics.HALF_SAT_MR_DESTROY_BACTERIA)
        ma_kills_bac = MacrophageDestroysBacterium(TBDynamics.RP_MA_DESTROY_BACTERIA, MACROPHAGE_ACTIVATED,
                                                   TBDynamics.HALF_SAT_MA_DESTROY_BACTERIA)
        events += [mr_kills_bac, ma_kills_bac]

        # bacterial translocation
        bed_translocation = CellTranslocationToLung(TBDynamics.RP_BACTERIA_BLOOD_TRANSLOCATION,
                                                    BACTERIUM_EXTRACELLULAR_DORMANT)
        events.append(bed_translocation)

        # Dendritic cell recruitment
        dc_recruit_lung = CellRecruitmentLung(TBDynamics.RP_DCI_RECRUITMENT, DENDRITIC_CELL_IMMATURE)
        events.append(dc_recruit_lung)

        # Dendritic cell death
        dci_death = CellDeath(TBDynamics.RP_DCI_DEATH, DENDRITIC_CELL_IMMATURE)
        dcm_death = CellDeath(TBDynamics.RP_DCM_DEATH, DENDRITIC_CELL_MATURE)
        events += [dci_death, dcm_death]

        # Dendritic Cell maturation
        dc_maturation = CellInfection(TBDynamics.RP_DCI_MATURATION, DENDRITIC_CELL_IMMATURE,
                                      TBDynamics.HALF_SAT_DCI_MATURATION)
        events.append(dc_maturation)

        # Dendritic cell translocation
        dc_translocation = CellTranslocationToLymph(TBDynamics.RP_DCM_TRANSLOCATION, DENDRITIC_CELL_MATURE)
        events.append(dc_translocation)

        # Macrophage recruitment
        mr_recruit_lung = CellRecruitmentLung(TBDynamics.RP_MR_RECRUIT_LUNG, MACROPHAGE_RESTING)
        mr_recruit_lymph = CellRecruitmentLymph(TBDynamics.RP_MR_RECRUIT_LYMPH, MACROPHAGE_RESTING)
        events += [mr_recruit_lung, mr_recruit_lymph]

        # Macrophage activation
        mr_activation_bac = CellActivation(TBDynamics.RP_MR_ACTIVATION_BY_BACTERIA,
                                           TBDynamics.HALF_SAT_MR_ACTIVATION_BY_BACTERIA, MACROPHAGE_RESTING,
                                           EXTRACELLULAR_BACTERIA)
        mr_activation_ta = CellActivation(TBDynamics.RP_MR_ACTIVATION_BY_TCELL,
                                          TBDynamics.HALF_SAT_MR_ACTIVATION_BY_TCELL,
                                          MACROPHAGE_RESTING, [T_CELL_ACTIVATED])
        events += [mr_activation_bac, mr_activation_ta]

        # Macrophage death
        mr_death = CellDeath(TBDynamics.RP_MR_DEATH, MACROPHAGE_RESTING)
        ma_death = CellDeath(TBDynamics.RP_MA_DEATH, MACROPHAGE_ACTIVATED)
        mi_death = CellDeath(TBDynamics.RP_MI_DEATH, MACROPHAGE_INFECTED)
        mi_burst = MacrophageBursting(TBDynamics.RP_MI_BURSTING, TBDynamics.SIGMOID_BIM_REPLICATION,
                                      TBDynamics.MACROPHAGE_CAPACITY)

        ta_kills_mi = TCellDestroysMacrophage(TBDynamics.RP_TA_KILL_MI, TBDynamics.HALF_SAT_TA_KILL_MI)
        events += [mr_death, ma_death, mi_death, mi_burst, ta_kills_mi]

        # Macrophage infection
        mac_infection = CellInfection(TBDynamics.RP_MR_INFECTION, MACROPHAGE_RESTING, TBDynamics.HALF_SAT_MR_INFECTION)
        events.append(mac_infection)

        # Macrophage translocation
        mi_translocation = CellTranslocationToLymph(TBDynamics.RP_MI_TRANSLOCATION, MACROPHAGE_INFECTED)
        events.append(mi_translocation)

        # T-cell recruitment
        tn_recruit = CellRecruitmentLymph(TBDynamics.RP_TN_RECRUIT, T_CELL_NAIVE)
        events.append(tn_recruit)

        # T-cell activation
        tn_activation = CellActivation(TBDynamics.RP_TCELL_ACTIVATION, TBDynamics.HALF_SAT_TCELL_ACTIVATION,
                                       T_CELL_NAIVE, [DENDRITIC_CELL_MATURE, MACROPHAGE_INFECTED])
        events.append(tn_activation)

        # T-cell translocation
        ta_translocation = TCellTranslocationToLungByInfection(TBDynamics.RP_TA_TRANSLOCATION)
        events.append(ta_translocation)

        # T-cell death
        tn_death = CellDeath(TBDynamics.RP_TN_DEATH, T_CELL_NAIVE)
        ta_death = CellDeath(TBDynamics.RP_TA_DEATH, T_CELL_ACTIVATED)
        events += [tn_death, ta_death]

        return events

    def _seed_network(self, params):
        # Seed attributes
        self._network.seed_pulmonary_attributes(params[PulmonaryNetwork.VENTILATION_SKEW],
                                                params[PulmonaryNetwork.PERFUSION_SKEW],
                                                params[PulmonaryNetwork.DRAINAGE_SKEW])
        # Seed initial compartments
        mac_recruit_lung = params[TBDynamics.RP_MR_RECRUIT_LUNG]
        mac_recruit_lymph = params[TBDynamics.RP_MR_RECRUIT_LYMPH]
        mac_death = params[TBDynamics.RP_MR_DEATH]
        dc_recruit = params[TBDynamics.RP_DCI_RECRUITMENT]
        dc_death = params[TBDynamics.RP_DCI_DEATH]
        tn_recruit = params[TBDynamics.RP_TN_RECRUIT]
        tn_death = params[TBDynamics.RP_TN_DEATH]
        self._network.seed_patches_by_rates({MACROPHAGE_RESTING: (mac_recruit_lung, mac_death),
                                             DENDRITIC_CELL_IMMATURE: (dc_recruit, dc_death)},
                                            {MACROPHAGE_RESTING: (mac_recruit_lymph, mac_death),
                                             T_CELL_NAIVE: (tn_recruit, tn_death)})

        # TODO - where to place bacteria (perfusion or fixed)?

    def _get_results(self):
        # TODO
        pass
