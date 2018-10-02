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
        self._ber_replication = Replication(BACTERIUM_EXTRACELLULAR_REPLICATING)
        self._bed_replication = Replication(BACTERIUM_EXTRACELLULAR_DORMANT)
        self._bim_replication = IntracellularBacterialReplication()
        events += [self._ber_replication, self._bed_replication, self._bim_replication]

        # Bacterial change
        self._bed_to_ber = BacteriumBecomesReplicating()
        self._ber_to_bed = BacteriumBecomesDormant()
        events += [self._bed_to_ber, self._ber_to_bed]

        # Bacterial destruction
        self._mr_kills_ber = MacrophageDestroysBacterium(MACROPHAGE_RESTING, BACTERIUM_EXTRACELLULAR_REPLICATING)
        self._mr_kills_bed = MacrophageDestroysBacterium(MACROPHAGE_RESTING, BACTERIUM_EXTRACELLULAR_DORMANT)
        self._ma_kills_ber = MacrophageDestroysBacterium(MACROPHAGE_ACTIVATED, BACTERIUM_EXTRACELLULAR_REPLICATING)
        self._ma_kills_bed = MacrophageDestroysBacterium(MACROPHAGE_ACTIVATED, BACTERIUM_EXTRACELLULAR_DORMANT)
        events += [self._mr_kills_ber, self._mr_kills_bed, self._ma_kills_ber, self._ma_kills_bed]

        # bacterial translocation
        self._bed_translocation = CellTranslocationToLung(BACTERIUM_EXTRACELLULAR_DORMANT)
        events.append(self._bed_translocation)

        # Dendritic cell recruitment
        self._dc_recruit_lung = CellRecruitmentLung(DENDRITIC_CELL_IMMATURE)
        events.append(self._dc_recruit_lung)

        # Dendritic cell death
        self._dci_death = CellDeath(DENDRITIC_CELL_IMMATURE)
        self._dcm_death = CellDeath(DENDRITIC_CELL_MATURE)
        events += [self._dci_death, self._dcm_death]

        # Dendritic Cell maturation
        self._dc_maturation = CellInfection(DENDRITIC_CELL_IMMATURE)
        events.append(self._dc_maturation)

        # Dendritic cell translocation
        self._dc_translocation = CellTranslocationToLymph(DENDRITIC_CELL_MATURE)
        events.append(self._dc_translocation)

        # Macrophage recruitment
        self._mr_recruit_lung = CellRecruitmentLung(MACROPHAGE_RESTING)
        self._mr_recruit_lymph = CellRecruitmentLymph(MACROPHAGE_RESTING)
        events += [self._mr_recruit_lung, self._mr_recruit_lymph]

        # Macrophage death
        self._mr_death = CellDeath(MACROPHAGE_RESTING)
        self._ma_death = CellDeath(MACROPHAGE_ACTIVATED)
        self._mi_death = CellDeath(MACROPHAGE_INFECTED)
        self._mi_burst = MacrophageBursting()
        self._ta_kills_mi = TCellDestroysMacrophage()
        events += [self._mr_death, self._ma_death, self._mi_death, self._mi_burst, self._ta_kills_mi]

        # Macrophage infection
        self._mac_infection = CellInfection(MACROPHAGE_RESTING)
        events.append(self._mac_infection)

        # Macrophage translocation
        self._mi_translocation = CellTranslocationToLymph(MACROPHAGE_INFECTED)
        events.append(self._mi_translocation)

        # Macrophage activation
        self._mr_activation_bac = CellActivation(MACROPHAGE_RESTING, EXTRACELLULAR_BACTERIA)
        self._mr_activation_ta = CellActivation(MACROPHAGE_RESTING, [T_CELL_ACTIVATED])
        events += [self._mr_activation_bac, self._mr_activation_ta]

        # T-cell recruitment
        self._tn_recruit = CellRecruitmentLymph(T_CELL_NAIVE)
        events.append(self._tn_recruit)

        # T-cell activation
        self._tn_activation = CellActivation(T_CELL_NAIVE, [DENDRITIC_CELL_MATURE, MACROPHAGE_INFECTED])
        events.append(self._tn_activation)

        # T-cell translocation
        self._ta_translocation = TCellTranslocationToLungByInfection()
        events.append(self._ta_translocation)

        # T-cell death
        self._tn_death = CellDeath(T_CELL_NAIVE)
        self._ta_death = CellDeath(T_CELL_ACTIVATED)
        events += [self._tn_death, self._ta_death]

        return events

    def _seed_events(self, params):
        self._ber_replication.set_reaction_parameter(params[TBDynamics.RP_REPLICATION_BER])
        self._bed_replication.set_reaction_parameter(params[TBDynamics.RP_REPLICATION_BED])
        self._bim_replication.set_parameters(params[TBDynamics.RP_REPLICATION_BIM],
                                             params[TBDynamics.SIGMOID_BIM_REPLICATION],
                                             params[TBDynamics.MACROPHAGE_CAPACITY])
        self._bed_to_ber.set_parameters(params[TBDynamics.RP_CHANGE_TO_REPLICATING],
                                        params[TBDynamics.SIGMOID_CHANGE_TO_REPLICATING],
                                        params[TBDynamics.HALF_SAT_CHANGE_TO_REPLICATING])
        self._ber_to_bed.set_parameters(params[TBDynamics.RP_CHANGE_TO_DORMANT],
                                        params[TBDynamics.SIGMOID_CHANGE_TO_DORMANT],
                                        params[TBDynamics.HALF_SAT_CHANGE_TO_DORMANT])
        for e in [self._mr_kills_bed, self._mr_kills_ber]:
            e.set_parameters(params[TBDynamics.RP_MR_DESTROY_BACTERIA], params[TBDynamics.HALF_SAT_MR_DESTROY_BACTERIA])
        for e in [self._ma_kills_bed, self._ma_kills_ber]:
            e.set_parameters(params[TBDynamics.RP_MA_DESTROY_BACTERIA], params[TBDynamics.HALF_SAT_MA_DESTROY_BACTERIA])
        self._bed_translocation.set_reaction_parameter(params[TBDynamics.RP_BACTERIA_BLOOD_TRANSLOCATION])

        self._dc_recruit_lung.set_reaction_parameter(params[TBDynamics.RP_DCI_RECRUITMENT])
        self._dci_death.set_reaction_parameter(params[TBDynamics.RP_DCI_DEATH])
        self._dcm_death.set_reaction_parameter(params[TBDynamics.RP_DCM_DEATH])
        self._dc_maturation.set_parameters(params[TBDynamics.RP_DCI_MATURATION],
                                           params[TBDynamics.HALF_SAT_DCI_MATURATION])
        self._dc_translocation.set_reaction_parameter(params[TBDynamics.RP_DCM_TRANSLOCATION])

        self._mr_recruit_lung.set_reaction_parameter(params[TBDynamics.RP_MR_RECRUIT_LUNG])
        self._mr_recruit_lymph.set_reaction_parameter(params[TBDynamics.RP_MR_RECRUIT_LYMPH])
        self._mr_activation_ta.set_parameters(params[TBDynamics.RP_MR_ACTIVATION_BY_TCELL],
                                              params[TBDynamics.HALF_SAT_MR_ACTIVATION_BY_TCELL])
        self._mr_activation_bac.set_parameters(params[TBDynamics.RP_MR_ACTIVATION_BY_BACTERIA],
                                               params[TBDynamics.HALF_SAT_MR_ACTIVATION_BY_BACTERIA])
        self._mr_death.set_reaction_parameter(params[TBDynamics.RP_MR_DEATH])
        self._ma_death.set_reaction_parameter(params[TBDynamics.RP_MA_DEATH])
        self._mi_death.set_reaction_parameter(params[TBDynamics.RP_MI_DEATH])
        self._mi_burst.set_parameters(params[TBDynamics.RP_MI_BURSTING],
                                      params[TBDynamics.SIGMOID_BIM_REPLICATION],
                                      params[TBDynamics.MACROPHAGE_CAPACITY])
        self._ta_kills_mi.set_parameters(params[TBDynamics.RP_TA_KILL_MI], params[TBDynamics.HALF_SAT_TA_KILL_MI])
        self._mac_infection.set_parameters(params[TBDynamics.RP_MR_INFECTION], params[TBDynamics.HALF_SAT_MR_INFECTION])
        self._mi_translocation.set_reaction_parameter(params[TBDynamics.RP_MI_TRANSLOCATION])

        self._tn_recruit.set_reaction_parameter(params[TBDynamics.RP_TN_RECRUIT])
        self._tn_activation.set_parameters(params[TBDynamics.RP_TCELL_ACTIVATION],
                                           params[TBDynamics.HALF_SAT_TCELL_ACTIVATION])
        self._ta_translocation.set_reaction_parameter(params[TBDynamics.RP_TA_TRANSLOCATION])
        self._tn_death.set_reaction_parameter(params[TBDynamics.RP_TN_DEATH])
        self._ta_death.set_reaction_parameter(params[TBDynamics.RP_TA_DEATH])

    def _seed_network(self, params):
        # Seed attributes
        self._network.seed_pulmonary_attributes(params)
        # Seed initial compartments
        mac_recruit_lung = params[TBDynamics.RP_MR_RECRUIT_LUNG]
        mac_recruit_lymph = params[TBDynamics.RP_MR_RECRUIT_LYMPH]
        mac_death = params[TBDynamics.RP_MR_DEATH]
        dc_recruit = params[TBDynamics.RP_DCI_RECRUITMENT]
        dc_death = params[TBDynamics.RP_DCI_DEATH]
        tn_recruit = params[TBDynamics.RP_TN_RECRUIT]
        tn_death = params[TBDynamics.RP_TN_DEATH]
        self._network.seed_alveolar_patches({MACROPHAGE_RESTING: (mac_recruit_lung, mac_death),
                                             DENDRITIC_CELL_IMMATURE: (dc_recruit, dc_death)})
        self._network.seed_lymph_patches({MACROPHAGE_RESTING: (mac_recruit_lymph, mac_death),
                                          T_CELL_NAIVE: (tn_recruit, tn_death)})

        # TODO - where to place bacteria (perfusion or fixed)?

    def _get_results(self):
        # TODO
        pass
