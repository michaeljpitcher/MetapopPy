from ..tbcompartments import *
from metapoppy.event import Event
import numpy


class CellInfection(Event):
    def __init__(self, cell_type):
        self._half_sat = 0
        self._cell_type = cell_type
        if cell_type == MACROPHAGE_RESTING:
            self._infected_cell_type = MACROPHAGE_INFECTED
            self._internalised_bac_type = BACTERIUM_INTRACELLULAR_MACROPHAGE
        elif cell_type == DENDRITIC_CELL_IMMATURE:
            self._infected_cell_type = DENDRITIC_CELL_MATURE
            self._internalised_bac_type = BACTERIUM_INTRACELLULAR_DENDRITIC
        else:
            raise Exception,"Invalid cell type"
        Event.__init__(self)

    def set_parameters(self, half_sat):
        self._half_sat = half_sat

    def _calculate_state_variable_at_patch(self, network, patch_id):
        total_bac = network.get_compartment_value(patch_id, BACTERIUM_EXTRACELLULAR_REPLICATING) + \
                    network.get_compartment_value(patch_id, BACTERIUM_EXTRACELLULAR_DORMANT)
        if not total_bac:
            return 0
        return network.get_compartment_value(patch_id, self._cell_type) * (float(total_bac) /
                                                                           (total_bac + self._half_sat))

    def perform(self, network, patch_id):
        replicating = network.get_compartment_value(patch_id, BACTERIUM_EXTRACELLULAR_REPLICATING)
        dormant = network.get_compartment_value(patch_id, BACTERIUM_EXTRACELLULAR_DORMANT)
        total_bacteria = replicating + dormant

        prob = numpy.array([replicating, dormant], dtype=numpy.dtype('float')) / total_bacteria
        bacteria_type_chosen = numpy.random.choice([BACTERIUM_EXTRACELLULAR_REPLICATING,
                                                    BACTERIUM_EXTRACELLULAR_DORMANT],
                                                   p=prob)

        changes = {self._cell_type: -1, self._infected_cell_type: 1, self._internalised_bac_type: 1,
                   bacteria_type_chosen: -1}
        network.update_patch(patch_id, changes)
