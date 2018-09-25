from metapoppy.dynamics import Dynamics
from .pulmonarynetwork import *
from tbmetapoppy.events import *
from tbcompartments import *


class TBDynamics(Dynamics):

    # Event Parameter keys
    MACROPHAGE_CAPACITY = 'macrophage_capacity'
    SIGMOID_BIM_REPLICATION = 'sigmoid_intracellular_bacteria_replication'
    RP_REPLICATION_BER = 'replicating_extracellular_bacteria_replication_rate'
    RP_REPLICATION_BED = 'dormant_extracellular_bacteria_replication_rate'
    RP_REPLICATION_BIM = 'intracellular_bacteria_replication_rate'
    RP_MR_DESTROY_BACTERIA = 'regular_macrophage_bacterial_destruction_rate'
    RP_MA_DESTROY_BACTERIA = 'activated_macrophage_bacterial_destruction_rate'
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
    # TODO - recruitment enhanced
    RP_DCM_TRANSLOCATION = 'dendrtic_cell_translocation_rate'

    RP_MACROPHAGE_ACTIVATION_BY_TCELL = 'macrophage_activation_by_tcell_rate'
    HALFSAT_MACROPHAGE_ACTIVATION_BY_TCELL = 'macrophage_activation_by_tcell_half_sat'
    RP_MACROPHAGE_ACTIVATION_BY_BACTERIA = 'macrophage_activation_by_bacteria_rate'
    HALFSAT_MACROPHAGE_ACTIVATION_BY_BACTERIA = 'macrophage_activation_by_bacteria_half_sat'
    RP_MACROPHAGE_INFECTION = 'macrophage_infection_rate'
    HALFSAT_MACROPHAGE_INFECTION = 'macrophage_infection_half_sat'
    RP_MR_DEATH = 'resting_macrophage_death_rate'
    RP_MI_DEATH = 'infected_macrophage_death_rate'
    RP_MA_DEATH = 'activated_macrophage_death_rate'
    RP_MACROPHAGE_BURSTING = 'macrophage_bursting_rate'
    RP_TA_KILL_MI = 't_cell_destroys_macrophage_rate'
    HALFSAT_TA_KILL_MI = 't_cell_destroys_macrophage_half_sat'
    RP_MACROPHAGE_RECRUITMENT_LUNG = 'macrophage_lung_recruitment_rate'
    RP_MACROPHAGE_RECRUITMENT_LYMPH = 'macrophage_lymph_recruitment_rate'
    # TODO - enhanced recruitment rates
    RP_MACROPHAGE_TRANSLOCATION = 'macrophage_translocation_rate'

    RP_TN_RECRUIT = 't_cell_recruitment_rate'
    # TODO - enhanced recruitment
    RP_TCELL_ACTIVATION = 't_cell_activation_rate'
    HALFSAT_TCELL_ACTIVATION = 't_cell_activation_half_sat'




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
        self._mi_translocation = CellTranslocationToLymph(DENDRITIC_CELL_MATURE)
        events.append(self._mi_translocation)

        # Macrophage activation
        self._mr_activation_bac = CellActivation(MACROPHAGE_RESTING, EXTRACELLULAR_BACTERIA)
        self._mr_activation_ta = CellActivation(MACROPHAGE_RESTING, [T_CELL_ACTIVATED])
        events += [self._mr_activation_bac, self._mr_activation_ta]

        # T-cell recruitment
        self._tn_recruit = CellRecruitmentLymph(T_CELL_NAIVE)
        events.append(self._tn_recruit)

        # T-cell activation
        self._tn_activation = CellActivation(T_CELL_ACTIVATED, [DENDRITIC_CELL_MATURE, MACROPHAGE_INFECTED])
        events.append(self._tn_activation)

        # T-cell death
        self._tn_death = CellDeath(T_CELL_NAIVE)
        self._ta_death = CellDeath(T_CELL_ACTIVATED)
        events += [self._tn_death, self._ta_death]

        # T-cell translocation
        self._ta_translocation = TCellTranslocationToLungByInfection()
        events.append(self._ta_translocation)

        return events

    def _seed_events(self, params):
        self._ber_replication.set_reaction_parameter(params[TBDynamics.RP_REPLICATION_BER])
        self._bed_replication.set_reaction_parameter(params[TBDynamics.RP_REPLICATION_BED])

    def _seed_network(self, params):
        # Seed attributes
        self._network.seed_pulmonary_attributes()
        # Seed initial compartments
        # TODO

    def _get_results(self):
        # TODO
        pass
