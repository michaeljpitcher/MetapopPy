from metapoppy_group_by_patch.network import TypedNetwork


class PulmonaryNetwork(TypedNetwork):

    # Patch types
    ALVEROLAR_PATCH = 'alveolar_patch'
    LYMPHATIC_PATCH = 'lymphatic_patch'

    # COMPARTMENTS
    BACTERIUM_EXTRACELLUAR_REPLICATING = 'bacterium_extracellular_replicating'
    BACTERIUM_EXTRACELLUAR_DORMANT = 'bacterium_extracellular_dormant'
    EXTRACELLULAR_BACTERIA = [BACTERIUM_EXTRACELLUAR_REPLICATING, BACTERIUM_EXTRACELLUAR_DORMANT]
    BACTERIUM_INTRACELLUAR_MACROPHAGE = 'bacterium_intracellular_macrophage'
    BACTERIUM_INTRACELLUAR_DENDRITIC = 'bacterium_intracellular_dendritic'
    INTRACELLULAR_BACTERIA = [BACTERIUM_INTRACELLUAR_DENDRITIC, BACTERIUM_INTRACELLUAR_MACROPHAGE]
    ALL_BACTERIA = INTRACELLULAR_BACTERIA + EXTRACELLULAR_BACTERIA

    MACROPHAGE_RESTING = 'macrophage_resting'
    MACROPHAGE_ACTIVATED = 'macrophage_activated'
    MACROPHAGE_INFECTED = 'macrophage_infected'
    ALL_MACROPHAGES = [MACROPHAGE_RESTING, MACROPHAGE_INFECTED, MACROPHAGE_ACTIVATED]

    DENDRITIC_CELL_IMMATURE = 'dendritic_cell_immature'
    DENDRITIC_CELL_MATURE = 'dendritic_cell_mature'
    ALL_DENDRITIC_CELLS = [DENDRITIC_CELL_IMMATURE, DENDRITIC_CELL_MATURE]

    T_CELL_NAIVE = 't_cell_naive'
    T_CELL_ACTIVATED = 't_cell_activated'
    ALL_T_CELLS = [T_CELL_NAIVE, T_CELL_ACTIVATED]

    ALL_COMPARTMENTS = ALL_BACTERIA + ALL_MACROPHAGES + ALL_DENDRITIC_CELLS + ALL_T_CELLS

    # ATTRIBUTES
    VENTILATION = 'ventilation'
    PERFUSION = 'perfusion'
    OXYGEN_TENSION = 'oxygen_tension'
    ALL_ALVEOLAR_ATTRIBUTES = [VENTILATION, PERFUSION, OXYGEN_TENSION]

    # Edge TODO - edge atts
    ALL_EDGE_ATTRIBUTES = []

    def __init__(self):
        TypedNetwork.__init__(self, PulmonaryNetwork.ALL_COMPARTMENTS, PulmonaryNetwork.ALL_ALVEOLAR_ATTRIBUTES,
                              PulmonaryNetwork.ALL_EDGE_ATTRIBUTES)

    def add_alveolar_patch(self, patch_id):
        self.add_node(patch_id)
        self.set_patch_type(patch_id, PulmonaryNetwork.ALVEROLAR_PATCH)

    def add_lymphatic_patch(self, patch_id):
        self.add_node(patch_id)
        self.set_patch_type(patch_id, PulmonaryNetwork.LYMPHATIC_PATCH)

    # TODO: build functions
