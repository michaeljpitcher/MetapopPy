from tbmetapoppy.tbpulmonarynetwork import TBPulmonaryNetwork
from metapoppy.event import PatchTypeEvent
from parameters import RATE, HALF_SAT
import numpy


class CellTranslocationToLymph(PatchTypeEvent):
    TRANSLOCATION_TO_LYMPH = '_translocation_to_lymph' + RATE

    def __init__(self, cell_type):
        self._cell_type = cell_type
        if cell_type in TBPulmonaryNetwork.INTERNAL_BACTERIA_FOR_CELL:
            self._internal_compartment = TBPulmonaryNetwork.INTERNAL_BACTERIA_FOR_CELL[cell_type]
        else:
            self._internal_compartment = None
        PatchTypeEvent.__init__(self, TBPulmonaryNetwork.ALVEOLAR_PATCH, [cell_type], [TBPulmonaryNetwork.DRAINAGE], [])

    def _define_parameter_keys(self):
        rp_key = self._cell_type + CellTranslocationToLymph.TRANSLOCATION_TO_LYMPH
        return rp_key, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._cell_type) * \
               network.get_attribute_value(patch_id, TBPulmonaryNetwork.DRAINAGE)

    def perform(self, network, patch_id):
        changes_lung = {self._cell_type: -1}
        changes_lymph = {self._cell_type: 1}
        # TODO - assumes only one lymph patch with PulmonaryNetwork.LYMPH_PATCH as patch ID
        lymph_id = TBPulmonaryNetwork.LYMPH_PATCH
        if self._internal_compartment:
            internals_to_transfer = int(round(float(
                                     network.get_compartment_value(patch_id, self._internal_compartment)) /
                                     network.get_compartment_value(patch_id, self._cell_type)))
            changes_lung[self._internal_compartment] = -1 * internals_to_transfer
            changes_lymph[self._internal_compartment] = internals_to_transfer
        network.update_patch(patch_id, changes_lung)
        network.update_patch(lymph_id, changes_lymph)


class CellTranslocationToLung(PatchTypeEvent):
    TRANSLOCATION_TO_LUNG = '_translocation_to_lung' + RATE

    def __init__(self, cell_type):
        self._cell_type = cell_type
        PatchTypeEvent.__init__(self, TBPulmonaryNetwork.LYMPH_PATCH, [cell_type], [], [])

    def _define_parameter_keys(self):
        rp_key = self._cell_type + CellTranslocationToLung.TRANSLOCATION_TO_LUNG
        return rp_key, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._cell_type)

    def perform(self, network, patch_id):
        edges = network.edges([patch_id],data=True)
        perfusions = [d[TBPulmonaryNetwork.PERFUSION] for _, _, d in edges]
        total_perfusion = sum(perfusions)

        lung_patch_id = numpy.random.choice([n for _,n,_ in edges], p=numpy.array(perfusions)/total_perfusion)
        network.update_patch(patch_id, {self._cell_type: -1})
        network.update_patch(lung_patch_id, {self._cell_type: 1})


class TCellTranslocationToLungByInfection(CellTranslocationToLung):
    def __init__(self):
        CellTranslocationToLung.__init__(self, TBPulmonaryNetwork.T_CELL_ACTIVATED)

    # TODO - this forces t-cells into the lung even if infection has been wiped out. State var maybe should depend
    #  on infection levels

    def perform(self, network, patch_id):
        infected_patches = network.infected_patches()
        if len(infected_patches) == 1:
            network.update_patch(patch_id, {self._cell_type: -1})
            network.update_patch(infected_patches[0], {self._cell_type: 1})
            return

        infections = [network.get_attribute_value(n, TBPulmonaryNetwork.PERFUSION) *
                      network.get_compartment_value(n, TBPulmonaryNetwork.EXTRACELLULAR_BACTERIA)
                      for n in infected_patches]
        total_infection = sum(infections)

        if not total_infection:
            return CellTranslocationToLung.perform(self, network, patch_id)

        lung_patch_id = numpy.random.choice(infected_patches, p=numpy.array(infections)/total_infection)
        network.update_patch(patch_id, {self._cell_type: -1})
        network.update_patch(lung_patch_id, {self._cell_type: 1})

# TODO translocation within lung
