from tbmetapoppy.tbpulmonarynetwork import TBPulmonaryNetwork
from metapoppy.event import PatchTypeEvent
from parameters import RATE, SIGMOID
import numpy


class Translocation(PatchTypeEvent):

    TRANSLOCATION_KEY = '_translocation_from_'

    def __init__(self, patch_type, cell_type):
        self._cell_type = cell_type
        PatchTypeEvent.__init__(self, patch_type, [cell_type], [], [])
        if cell_type in TBPulmonaryNetwork.INTERNAL_BACTERIA_FOR_CELL:
            self._internal_compartment = TBPulmonaryNetwork.INTERNAL_BACTERIA_FOR_CELL[cell_type]
        else:
            self._internal_compartment = None

    def _define_parameter_keys(self):
        return self._cell_type + Translocation.TRANSLOCATION_KEY + self._patch_type + RATE, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._cell_type)

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


class TranslocationLungToLymph(Translocation):
    def __init__(self, cell_type):
        Translocation.__init__(self, TBPulmonaryNetwork.ALVEOLAR_PATCH, cell_type)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._cell_type) * \
               network.get_attribute_value(patch_id, TBPulmonaryNetwork.DRAINAGE)


class TranslocationLymphToLung(Translocation):
    def __init__(self, cell_type):
        Translocation.__init__(self, TBPulmonaryNetwork.LYMPH_PATCH, cell_type)

    def _define_parameter_keys(self):
        self._sigmoid_key = self._cell_type + Translocation.TRANSLOCATION_KEY + self._patch_type + SIGMOID
        return self._cell_type + Translocation.TRANSLOCATION_KEY + self._patch_type + RATE, [self._sigmoid_key]

    def _calculate_state_variable_at_patch(self, network, patch_id):
        cells = network.get_compartment_value(patch_id, self._cell_type)
        if not cells:
            return 0
        infected_patches = network.infected_patches()
        cytokine_count_lung = sum([network.get_edge_data(patch_id, n)[TBPulmonaryNetwork.CYTOKINE] for n in infected_patches])
        cytokine_count_lymph = network.get_compartment_value(patch_id, TBPulmonaryNetwork.MACROPHAGE_INFECTED)
        # Catch to avoid / 0 errors
        if not cytokine_count_lymph and not cytokine_count_lung:
            return 0
        return network.get_compartment_value(patch_id, TBPulmonaryNetwork.T_CELL_ACTIVATED) * \
               (float(cytokine_count_lung) ** self._parameters[self._sigmoid_key] /
                (cytokine_count_lung ** self._parameters[self._sigmoid_key] +
                 cytokine_count_lymph ** self._parameters[self._sigmoid_key]))

    def _choose_neighbour(self, network, patch_id):
        # Choosing based on infection and perfusion
        infected_patches = network.infected_patches()
        edges = {n: network.get_edge_data(patch_id, n) for n in infected_patches}
        neighbours = []
        vals = []
        total = 0
        for k,v in edges.iteritems():
            neighbours.append(k)
            val = v[TBPulmonaryNetwork.CYTOKINE] * v[TBPulmonaryNetwork.PERFUSION]
            vals.append(val)
            total += val
        # Choose a neighbour based on the values
        return neighbours[numpy.random.choice(range(len(neighbours)), p=numpy.array(vals)/total)]

