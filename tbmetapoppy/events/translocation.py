from tbmetapoppy.tbpulmonaryenvironment import TBPulmonaryEnvironment
from metapoppy.event import PatchTypeEvent
from parameters import RATE, SIGMOID, HALF_SAT
import numpy


class TranslocationLungToLymph(PatchTypeEvent):
    TRANSLOCATION_KEY = '_translocation_from_'

    def __init__(self, cell_type):
        self._cell_type = cell_type
        dep_comps = [cell_type]

        if cell_type in TBPulmonaryEnvironment.INTERNAL_BACTERIA_FOR_CELL:
            self._internal_compartment = TBPulmonaryEnvironment.INTERNAL_BACTERIA_FOR_CELL[cell_type]
        else:
            self._internal_compartment = None

        PatchTypeEvent.__init__(self, TBPulmonaryEnvironment.ALVEOLAR_PATCH, dep_comps,
                                [TBPulmonaryEnvironment.DRAINAGE], [])

    def _define_parameter_keys(self):
        return self._cell_type + TranslocationLungToLymph.TRANSLOCATION_KEY + self._patch_type + RATE, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._cell_type) * \
               network.get_attribute_value(patch_id, TBPulmonaryEnvironment.DRAINAGE)

    def perform(self, network, patch_id):
        neighbour = TBPulmonaryEnvironment.LYMPH_PATCH
        changes_from = {self._cell_type: -1}
        changes_to = {self._cell_type: 1}
        if self._internal_compartment:
            internals_to_transfer = int(round(float(
                network.get_compartment_value(patch_id, self._internal_compartment)) /
                network.get_compartment_value(patch_id, self._cell_type)))
            changes_from[self._internal_compartment] = -1 * internals_to_transfer
            changes_to[self._internal_compartment] = internals_to_transfer
        network.update_patch(patch_id, changes_from)
        network.update_patch(neighbour, changes_to)


class TranslocationLymphToLungCytokine(PatchTypeEvent):
    TRANSLOCATION_KEY = '_translocation_from_'

    def __init__(self, cell_type):
        self._cell_type = cell_type
        dep_comps = [cell_type, TBPulmonaryEnvironment.MACROPHAGE_INFECTED]

        PatchTypeEvent.__init__(self, TBPulmonaryEnvironment.LYMPH_PATCH, dep_comps,
                                [], [TBPulmonaryEnvironment.CYTOKINE])

    def _define_parameter_keys(self):
        self._sigmoid_key = self._cell_type + TranslocationLungToLymph.TRANSLOCATION_KEY + self._patch_type + '_by_' + \
               TBPulmonaryEnvironment.CYTOKINE + SIGMOID
        return self._cell_type + TranslocationLungToLymph.TRANSLOCATION_KEY + self._patch_type + '_by_' + \
               TBPulmonaryEnvironment.CYTOKINE + RATE, [self._sigmoid_key]

    def _calculate_state_variable_at_patch(self, network, patch_id):
        cells = network.get_compartment_value(patch_id, self._cell_type)
        if not cells:
            return 0
        infected_patches = network.infected_patches()
        cytokine_count_lung = sum([network.get_edge_data(patch_id, n)[TBPulmonaryEnvironment.CYTOKINE] for n in
                                   infected_patches])
        cytokine_count_lymph = network.get_compartment_value(patch_id, TBPulmonaryEnvironment.MACROPHAGE_INFECTED)
        # Catch to avoid / 0 errors
        if not cytokine_count_lymph and not cytokine_count_lung:
            return 0
        return cells * (float(cytokine_count_lung) ** self._parameters[self._sigmoid_key] /
                        (cytokine_count_lung ** self._parameters[self._sigmoid_key] +
                         cytokine_count_lymph ** self._parameters[self._sigmoid_key]))

    def perform(self, network, patch_id):
        infected_patches = network.infected_patches()

        edges = {n: network[patch_id][n][TBPulmonaryEnvironment.PERFUSION] for n in infected_patches}
        neighbours = []
        vals = []
        total = 0
        for k, v in edges.iteritems():
            neighbours.append(k)
            val = network.get_compartment_value(k, TBPulmonaryEnvironment.MACROPHAGE_INFECTED) * v
            vals.append(val)
            total += val

        # Choose a neighbour based on the values
        neighbour = neighbours[numpy.random.choice(range(len(neighbours)), p=numpy.array(vals) / total)]

        network.update_patch(patch_id, {self._cell_type: -1})
        network.update_patch(neighbour, {self._cell_type: 1})


class TranslocationLymphToLungDendritic(PatchTypeEvent):
    TRANSLOCATION_KEY = '_translocation_from_'

    def __init__(self, cell_type):
        self._cell_type = cell_type
        dep_comps = [cell_type, TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE]

        PatchTypeEvent.__init__(self, TBPulmonaryEnvironment.LYMPH_PATCH, dep_comps,
                                [], [])

    def _define_parameter_keys(self):
        self._sigmoid_key = self._cell_type + TranslocationLungToLymph.TRANSLOCATION_KEY + self._patch_type + '_by_' + \
                            TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE + SIGMOID
        self._half_sat_key = self._cell_type + TranslocationLungToLymph.TRANSLOCATION_KEY + self._patch_type + '_by_' + \
                            TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE + HALF_SAT
        return self._cell_type + TranslocationLungToLymph.TRANSLOCATION_KEY + self._patch_type + '_by_' + \
               TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE + RATE, [self._sigmoid_key, self._half_sat_key]

    def _calculate_state_variable_at_patch(self, network, patch_id):
        cells = network.get_compartment_value(patch_id, self._cell_type)
        if not cells:
            return 0
        dm = network.get_compartment_value(patch_id, TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE)
        # Catch to avoid / 0 errors
        if not dm:
            return 0
        sig = self._parameters[self._sigmoid_key]
        return cells * (float(dm)**sig / (dm**sig + self._parameters[self._half_sat_key]**sig))

    def perform(self, network, patch_id):
        infected_patches = network.infected_patches()

        edges = {n: network[patch_id][n][TBPulmonaryEnvironment.PERFUSION] for n in infected_patches}
        neighbours = []
        vals = []
        total = 0
        for k, v in edges.iteritems():
            neighbours.append(k)
            val = network.get_compartment_value(k, TBPulmonaryEnvironment.MACROPHAGE_INFECTED) * v
            vals.append(val)
            total += val

        if total == 0:
            return

        # Choose a neighbour based on the values
        neighbour = neighbours[numpy.random.choice(range(len(neighbours)), p=numpy.array(vals) / total)]

        network.update_patch(patch_id, {self._cell_type: -1})
        network.update_patch(neighbour, {self._cell_type: 1})


class TranslocationLymphToLungBlood(PatchTypeEvent):
    TRANSLOCATION_KEY = '_translocation_from_'
    BLOOD = 'blood'

    def __init__(self, cell_type):
        self._cell_type = cell_type
        dep_comps = [cell_type, TBPulmonaryEnvironment.SOLID_CASEUM]

        PatchTypeEvent.__init__(self, TBPulmonaryEnvironment.LYMPH_PATCH, dep_comps,
                                [], [])

    def _define_parameter_keys(self):
        self._half_sat_key = self._cell_type + TranslocationLungToLymph.TRANSLOCATION_KEY + self._patch_type + '_by_' + \
                            TranslocationLymphToLungBlood.BLOOD + HALF_SAT
        return self._cell_type + TranslocationLungToLymph.TRANSLOCATION_KEY + self._patch_type + '_by_' + \
                            TranslocationLymphToLungBlood.BLOOD + RATE, [self._half_sat_key]

    def _calculate_state_variable_at_patch(self, network, patch_id):
        cells = network.get_compartment_value(patch_id, self._cell_type)
        if not cells:
            return 0
        cas = network.get_compartment_value(patch_id, TBPulmonaryEnvironment.SOLID_CASEUM)
        return cells * (1 - (float(cas) / (cas + self._parameters[self._half_sat_key])))

    def perform(self, network, patch_id):
        # TODO - this works on the basis of sum of all perfusion == 1, needs amending if perfusion changes
        r = numpy.random.random()
        edges = network[patch_id]
        total = 0
        for k, v in edges.iteritems():
            total += v[TBPulmonaryEnvironment.PERFUSION]
            if total >= r:
                network.update_patch(patch_id, {self._cell_type: -1})
                network.update_patch(k, {self._cell_type: 1})
                return
