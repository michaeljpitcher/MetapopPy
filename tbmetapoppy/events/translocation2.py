from tbmetapoppy.tbpulmonaryenvironment import TBPulmonaryEnvironment
from metapoppy.event import PatchTypeEvent
from parameters import RATE, SIGMOID, HALF_SAT
import numpy


class Translocation(PatchTypeEvent):

    TRANSLOCATION_KEY = '_translocation_from_'

    def __init__(self, patch_type, cell_type, influencing_comps=None):
        self._moving_compartment = cell_type
        dep_comps = [cell_type]
        if influencing_comps:
            dep_comps += influencing_comps
        PatchTypeEvent.__init__(self, patch_type, dep_comps, [], [])
        if cell_type in TBPulmonaryEnvironment.INTERNAL_BACTERIA_FOR_CELL:
            self._internal_compartment = TBPulmonaryEnvironment.INTERNAL_BACTERIA_FOR_CELL[cell_type]
        else:
            self._internal_compartment = None

    def _define_parameter_keys(self):
        return self._moving_compartment + Translocation.TRANSLOCATION_KEY + self._patch_type + RATE, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._moving_compartment)

    def perform(self, network, patch_id):
        neighbour = self._choose_neighbour(network, patch_id)
        self._translocate(network, patch_id, neighbour)

    def _choose_neighbour(self, network, patch_id):
        edges = network[patch_id]
        if len(edges) == 1:
            return edges.keys()[0]
        else:
            return numpy.random.choice(edges.keys())

    def _translocate(self, network, patch_id, neighbour):
        changes_from = {self._moving_compartment: -1}
        changes_to = {self._moving_compartment: 1}
        if self._internal_compartment:
            internals_to_transfer = int(round(float(
                network.get_compartment_value(patch_id, self._internal_compartment)) /
                                              network.get_compartment_value(patch_id, self._moving_compartment)))
            changes_from[self._internal_compartment] = -1 * internals_to_transfer
            changes_to[self._internal_compartment] = internals_to_transfer
        network.update_patch(patch_id, changes_from)
        network.update_patch(neighbour, changes_to)


class TranslocationLungToLymph(Translocation):
    def __init__(self, cell_type):
        Translocation.__init__(self, TBPulmonaryEnvironment.ALVEOLAR_PATCH, cell_type)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._moving_compartment) * \
               network.get_attribute_value(patch_id, TBPulmonaryEnvironment.DRAINAGE)


class TranslocationLymphToLung(Translocation):
    TRANSLOCATION_CASEUM_HALF_SAT = 'caseum_half_sat'

    def __init__(self, cell_type):
        Translocation.__init__(self, TBPulmonaryEnvironment.LYMPH_PATCH, cell_type, [TBPulmonaryEnvironment.SOLID_CASEUM])

    def _define_parameter_keys(self):
        self._caseum_key = TranslocationLymphToLung.TRANSLOCATION_CASEUM_HALF_SAT
        return self._moving_compartment + Translocation.TRANSLOCATION_KEY + self._patch_type + RATE, [self._caseum_key]

    def _calculate_state_variable_at_patch(self, network, patch_id):
        caseum = network.get_compartment_value(patch_id, TBPulmonaryEnvironment.SOLID_CASEUM)
        return network.get_compartment_value(patch_id, self._moving_compartment) * \
               (1 - (float(caseum)/(caseum+self._parameters[self._caseum_key])))

    def _choose_neighbour(self, network, patch_id):
        # Choosing based on perfusion
        edges = network[patch_id]
        neighbours = []
        vals = []
        # total = 0
        for k, v in edges.iteritems():
            neighbours.append(k)
            val = v[TBPulmonaryEnvironment.PERFUSION]
            vals.append(float(val))
            # total += val
        # Choose a neighbour
        n = numpy.random.choice(neighbours, p=vals)
        return n


class TCellTranslocationLymphToLung(Translocation):
    def __init__(self, cell_type):
        Translocation.__init__(self, TBPulmonaryEnvironment.LYMPH_PATCH, cell_type,
                               [TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE,
                                TBPulmonaryEnvironment.MACROPHAGE_INFECTED])

    def _define_parameter_keys(self):
        self._sigmoid_key = self._moving_compartment + Translocation.TRANSLOCATION_KEY + self._patch_type + SIGMOID
        self._half_sat_key = self._moving_compartment + Translocation.TRANSLOCATION_KEY + self._patch_type + HALF_SAT
        return self._moving_compartment + Translocation.TRANSLOCATION_KEY + self._patch_type + RATE, \
               [self._sigmoid_key, self._half_sat_key]

    def _calculate_state_variable_at_patch(self, network, patch_id):
        cells = network.get_compartment_value(patch_id, self._moving_compartment)
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

    # def _calculate_state_variable_at_patch(self, network, patch_id):
    #     cells = network.get_compartment_value(patch_id, self._moving_compartment)
    #     if not cells:
    #         return 0
    #     dm = network.get_compartment_value(patch_id, TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE)
    #     # Catch to avoid / 0 errors
    #     if not dm:
    #         return 0
    #     sig = self._parameters[self._sigmoid_key]
    #     return cells * (float(dm)**sig / (dm**sig + self._parameters[self._half_sat_key]**sig))

    def _choose_neighbour(self, network, patch_id):
        # Choosing based on infection and perfusion
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
        return neighbours[numpy.random.choice(range(len(neighbours)), p=numpy.array(vals)/total)]
