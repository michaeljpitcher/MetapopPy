from tbmetapoppy import *
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

        lymph_id = next(network.edges([patch_id], data=True))[1]

        if self._cell_type == TBDynamics.MACROPHAGE_INFECTED:
            bac_to_transfer = int(round(float(
                network.get_compartment_value(patch_id, TBDynamics.BACTERIUM_INTRACELLULAR_MACROPHAGE) /
                network.get_compartment_value(patch_id, TBDynamics.MACROPHAGE_INFECTED))))
            changes_lung[TBDynamics.BACTERIUM_INTRACELLULAR_MACROPHAGE] = -1 * bac_to_transfer
            changes_lymph[TBDynamics.BACTERIUM_INTRACELLULAR_MACROPHAGE] = bac_to_transfer
        elif self._cell_type == TBDynamics.DENDRITIC_CELL_MATURE:
            changes_lung[TBDynamics.BACTERIUM_INTRACELLULAR_DENDRITIC] = -1
            changes_lymph[TBDynamics.BACTERIUM_INTRACELLULAR_DENDRITIC] = 1

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
        total_perfusion = sum(d[PulmonaryNetwork.PERFUSION] for _,_,d in edges)

        lung_patch_id = numpy.random.choice(edges,
                                p=numpy.array([d[PulmonaryNetwork.PERFUSION] for _, _, d in edges])/total_perfusion)[1]
        network.update_patch(patch_id, {self._cell_type: -1})
        network.update_patch(lung_patch_id, {self._cell_type: 1})

# TODO translocation within lung
