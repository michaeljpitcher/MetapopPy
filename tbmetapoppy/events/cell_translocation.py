from tbmetapoppy.pulmonarynetwork import PulmonaryNetwork
from ..tbcompartments import *
from metapoppy.event import PatchTypeEvent
import numpy


class CellTranslocationToLymph(PatchTypeEvent):
    def __init__(self, cell_type):
        self._cell_type = cell_type
        PatchTypeEvent.__init__(self, PulmonaryNetwork.ALVEOLAR_PATCH)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._cell_type) * \
               network.get_attribute_value(patch_id, PulmonaryNetwork.DRAINAGE)

    def perform(self, network, patch_id):
        changes_lung = {self._cell_type: -1}
        changes_lymph = {self._cell_type: 1}

        # TODO - assumes only one lymph patch with PulmonaryNetwork.LYMPH_PATCH as patch ID
        lymph_id = PulmonaryNetwork.LYMPH_PATCH

        if self._cell_type == MACROPHAGE_INFECTED:
            bac_to_transfer = int(round(float(
                network.get_compartment_value(patch_id, BACTERIUM_INTRACELLULAR_MACROPHAGE)) /
                network.get_compartment_value(patch_id, MACROPHAGE_INFECTED)))
            changes_lung[BACTERIUM_INTRACELLULAR_MACROPHAGE] = -1 * bac_to_transfer
            changes_lymph[BACTERIUM_INTRACELLULAR_MACROPHAGE] = bac_to_transfer
        elif self._cell_type == DENDRITIC_CELL_MATURE:
            changes_lung[BACTERIUM_INTRACELLULAR_DENDRITIC] = -1
            changes_lymph[BACTERIUM_INTRACELLULAR_DENDRITIC] = 1

        network.update_patch(patch_id, changes_lung)
        network.update_patch(lymph_id, changes_lymph)


class CellTranslocationToLung(PatchTypeEvent):
    def __init__(self, cell_type):
        self._cell_type = cell_type
        PatchTypeEvent.__init__(self, PulmonaryNetwork.LYMPH_PATCH)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._cell_type)

    def perform(self, network, patch_id):
        edges = network.edges([patch_id],data=True)
        perfusions = [d[PulmonaryNetwork.PERFUSION] for _, _, d in edges]
        total_perfusion = sum(perfusions)

        lung_patch_id = numpy.random.choice([n for _,n,_ in edges], p=numpy.array(perfusions)/total_perfusion)
        network.update_patch(patch_id, {self._cell_type: -1})
        network.update_patch(lung_patch_id, {self._cell_type: 1})


class TCellTranslocationToLungByInfection(CellTranslocationToLung):
    def __init__(self):
        CellTranslocationToLung.__init__(self, T_CELL_ACTIVATED)

    def perform(self, network, patch_id):
        edges = network.edges([patch_id],data=True)

        infections = [d[PulmonaryNetwork.PERFUSION] *
                      (network.get_compartment_value(n, BACTERIUM_EXTRACELLULAR_DORMANT) +
                       network.get_compartment_value(n, BACTERIUM_EXTRACELLULAR_REPLICATING)) for _, n, d in edges]

        total_infection = sum(infections)

        if not total_infection:
            return CellTranslocationToLung.perform(self, network, patch_id)

        lung_patch_id = numpy.random.choice([n for _,n,_ in edges], p=numpy.array(infections)/total_infection)
        network.update_patch(patch_id, {self._cell_type: -1})
        network.update_patch(lung_patch_id, {self._cell_type: 1})

# TODO translocation within lung
