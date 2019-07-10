from metapoppy.dynamics import Dynamics
from tbmetapoppy.events import *


class TBDynamics(Dynamics):
    # Attribute seeding params

    IC_BAC_LOCATION = 'initial_bacterial_load_location'
    RANDOM = 'random'
    IC_BER_LOAD = 'initial_bacterial_load_replicating'
    IC_BED_LOAD = 'initial_bacterial_load_dormant'

    def __init__(self, network_config):
        self._lung_recruit_keys = {}
        self._lymph_recruit_keys = {}
        self._death_rates = {}

        self._perf_seed = {}

        # Build network
        pulmonary_network = TBPulmonaryEnvironment(network_config)
        Dynamics.__init__(self, pulmonary_network)

    def output_positions(self, filename):
        self._prototype_network.output_positions(filename)

    def _create_events(self):
        """
        Create all TB related events, and assign relevant parameters to them
        :return:
        """
        events = []

        # Bacteria replication
        ber_replication = Replication(TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_REPLICATING)
        bed_replication = Replication(TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT)
        bim_replication = IntracellularBacterialReplication()

        events += [ber_replication, bed_replication, bim_replication]

        # Bacterial change
        bed_to_ber = BacteriumChangeStateThroughOxygen(False)
        ber_to_bed = BacteriumChangeStateThroughOxygen(True)
        events += [bed_to_ber, ber_to_bed]

        # Bacterial translocation - lymph to lung
        bed_translocation = BacteriaTranslocationLymphToLung(TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT)
        events.append(bed_translocation)

        # Dendritic cell recruitment
        dc_recruit_lung = StandardCellRecruitmentLung(TBPulmonaryEnvironment.DENDRITIC_CELL_IMMATURE)
        self._lung_recruit_keys[TBPulmonaryEnvironment.DENDRITIC_CELL_IMMATURE] = dc_recruit_lung.reaction_parameter()
        dc_recruit_lung_enh = EnhancedCellRecruitmentLung(TBPulmonaryEnvironment.DENDRITIC_CELL_IMMATURE)
        events += [dc_recruit_lung, dc_recruit_lung_enh]

        # Dendritic cell death
        dci_death = CellDeath(TBPulmonaryEnvironment.DENDRITIC_CELL_IMMATURE)
        self._death_rates[TBPulmonaryEnvironment.DENDRITIC_CELL_IMMATURE] = dci_death.reaction_parameter()
        dcm_death = InfectedCellDeath(TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE)
        events += [dci_death, dcm_death]

        # Dendritic Cell maturation
        dc_maturation = CellIngestBacterium(TBPulmonaryEnvironment.DENDRITIC_CELL_IMMATURE)
        events.append(dc_maturation)

        # Dendritic cell translocation
        dc_translocation = TranslocationLungToLymph(TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE)
        events.append(dc_translocation)

        # Macrophage recruitment
        mr_recruit_lung = StandardCellRecruitmentLung(TBPulmonaryEnvironment.MACROPHAGE_RESTING)
        self._lung_recruit_keys[TBPulmonaryEnvironment.MACROPHAGE_RESTING] = mr_recruit_lung.reaction_parameter()
        mr_recruit_lymph = StandardCellRecruitmentLymph(TBPulmonaryEnvironment.MACROPHAGE_RESTING)
        self._lymph_recruit_keys[TBPulmonaryEnvironment.MACROPHAGE_RESTING] = mr_recruit_lymph.reaction_parameter()
        mr_recruit_lung_enhanced = EnhancedCellRecruitmentLung(TBPulmonaryEnvironment.MACROPHAGE_RESTING)
        mr_recruit_lymph_enhanced = EnhancedCellRecruitmentLymph(TBPulmonaryEnvironment.MACROPHAGE_RESTING)
        
        events += [mr_recruit_lung, mr_recruit_lymph, mr_recruit_lung_enhanced, mr_recruit_lymph_enhanced]

        # Macrophage activation
        mr_activation_bac = CellActivation(TBPulmonaryEnvironment.MACROPHAGE_RESTING,
                                           TBPulmonaryEnvironment.EXTRACELLULAR_BACTERIA)
        mr_activation_ta = CellActivation(TBPulmonaryEnvironment.MACROPHAGE_RESTING,
                                          [TBPulmonaryEnvironment.T_CELL_ACTIVATED])
        events += [mr_activation_bac, mr_activation_ta]

        # Macrophage death
        mr_death = CellDeath(TBPulmonaryEnvironment.MACROPHAGE_RESTING)
        self._death_rates[TBPulmonaryEnvironment.MACROPHAGE_RESTING] = mr_death.reaction_parameter()
        ma_death = CellDeath(TBPulmonaryEnvironment.MACROPHAGE_ACTIVATED)
        mi_death = InfectedCellDeath(TBPulmonaryEnvironment.MACROPHAGE_INFECTED)
        mi_burst = MacrophageBursting()

        ta_kills_mi = TCellDestroysMacrophage()
        events += [mr_death, ma_death, mi_death, mi_burst, ta_kills_mi]

        # Macrophage ingest bacterium
        mr_ingest_bac = CellIngestBacterium(TBPulmonaryEnvironment.MACROPHAGE_RESTING)
        ma_ingest_bac = CellIngestBacterium(TBPulmonaryEnvironment.MACROPHAGE_ACTIVATED)
        events += [mr_ingest_bac, ma_ingest_bac]

        # Macrophage translocation
        mi_translocation = TranslocationLungToLymph(TBPulmonaryEnvironment.MACROPHAGE_INFECTED)
        events.append(mi_translocation)

        # T-cell recruitment
        tn_recruit = StandardCellRecruitmentLymph(TBPulmonaryEnvironment.T_CELL_NAIVE)
        self._lymph_recruit_keys[TBPulmonaryEnvironment.T_CELL_NAIVE] = tn_recruit.reaction_parameter()
        events.append(tn_recruit)

        # T-cell enhanced recruitment
        tn_recruit_enhanced = EnhancedTCellRecruitmentLymph()
        events.append(tn_recruit_enhanced)

        # T-cell activation
        tn_activation = CellActivation(TBPulmonaryEnvironment.T_CELL_NAIVE,
                                       [TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE,
                                        TBPulmonaryEnvironment.MACROPHAGE_INFECTED])
        events.append(tn_activation)

        ta_rep = Replication(TBPulmonaryEnvironment.T_CELL_ACTIVATED)
        events.append(ta_rep)

        # T-cell translocation
        # ta_translocation = TCellTranslocationToLungByInfection()
        ta_translocation = TCellTranslocationLymphToLung(TBPulmonaryEnvironment.T_CELL_ACTIVATED)
        events.append(ta_translocation)

        # T-cell death
        tn_death = CellDeath(TBPulmonaryEnvironment.T_CELL_NAIVE)
        self._death_rates[TBPulmonaryEnvironment.T_CELL_NAIVE] = tn_death.reaction_parameter()
        ta_death = CellDeath(TBPulmonaryEnvironment.T_CELL_ACTIVATED)
        events += [tn_death, ta_death]

        return events

    def _get_initial_patch_seeding(self, params):
        """
        Get the initial values for patches. Only needs bacteria - other cells will be added when patches become
        active
        :param params:
        :return:
        """
        # Calculate the initial patch attribute values
        patch_seeding = self._network.calculate_pulmonary_attribute_values(params)

        # Bacteria
        # Ventilation based - assumes sum of ventilation values = 1.0
        initial_bac_patch = None

        if params[TBDynamics.IC_BAC_LOCATION] == TBDynamics.RANDOM:
            r = numpy.random.random()
            count = 0
            for p in self._network.get_patches_by_type(TBPulmonaryEnvironment.ALVEOLAR_PATCH):
                count += patch_seeding[p][TypedEnvironment.ATTRIBUTES][TBPulmonaryEnvironment.VENTILATION]
                if count >= r:
                    initial_bac_patch = p
                    break
        else:
            initial_bac_patch = params[TBDynamics.IC_BAC_LOCATION]


        patch_seeding[initial_bac_patch][TypedEnvironment.COMPARTMENTS] = \
            {TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_REPLICATING: params[TBDynamics.IC_BER_LOAD],
             TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT: params[TBDynamics.IC_BED_LOAD]}

        return patch_seeding

    def _patch_is_active(self, patch_id):
        """
        Determine if the given patch is active (from the network). Alveolar patches only become active when they contain
        bacteria.
        :param patch_id:
        :return:
        """
        return patch_id == TBPulmonaryEnvironment.LYMPH_PATCH or \
                sum([self._network.get_compartment_value(patch_id, n) for n in TBPulmonaryEnvironment.BACTERIA]) > 0

    def _seed_activated_patch(self, patch_id, params):
        """
        When a patch becomes activated, seed it with the equilibrium level of immune cells
        :param patch_id:
        :param params:
        :return:
        """

        death_rates = {c: params[self._death_rates[c]] for c in
                       [TBPulmonaryEnvironment.MACROPHAGE_RESTING, TBPulmonaryEnvironment.T_CELL_NAIVE,
                        TBPulmonaryEnvironment.DENDRITIC_CELL_IMMATURE]}

        if patch_id == TBPulmonaryEnvironment.LYMPH_PATCH:
            # Attribute seeding
            patch_seeding = {}

            lymph_recruit_rates = {c: params[self._lymph_recruit_keys[c]] for c in
                                   [TBPulmonaryEnvironment.MACROPHAGE_RESTING, TBPulmonaryEnvironment.T_CELL_NAIVE]}
            patch_seeding[TypedEnvironment.COMPARTMENTS] = \
                {c: int(round(lymph_recruit_rates[c] / death_rates[c])) for c in lymph_recruit_rates}

        else:
            # Attribute seeding
            patch_seeding = {}

            lung_recruit_rates = {c: params[self._lung_recruit_keys[c]] for c in
                                  [TBPulmonaryEnvironment.MACROPHAGE_RESTING,
                                   TBPulmonaryEnvironment.DENDRITIC_CELL_IMMATURE]}

            perf = self._network.get_attribute_value(patch_id, TBPulmonaryEnvironment.PERFUSION)

            patch_seeding[TypedEnvironment.COMPARTMENTS] = \
                {c: int(round(perf * (lung_recruit_rates[c] / death_rates[c]))) for c in lung_recruit_rates}

        return patch_seeding

    def _build_network(self, params):
        pass

    def _get_initial_edge_seeding(self, params):
        pass

    def _end_simulation(self, t):
        # DEBUG - remove
        return any([i for i in self._network._infected_patches if
                    self._network.get_compartment_value(i, TBPulmonaryEnvironment.BACTERIA) >= 20000])
