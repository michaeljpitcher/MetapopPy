from ..tbcompartments import *
from metapoppy.event import Event
import numpy


class CellInfection(Event):
    def __init__(self, infection_rate_key, cell_type, half_sat_key):
        self._half_sat_key = half_sat_key
        self._cell_type = cell_type
        self._infected_cell_type = INFECTED_CELL[self._cell_type]
        self._internalised_bac_type = INTERNAL_BACTERIA_FOR_CELL[self._infected_cell_type]
        Event.__init__(self, infection_rate_key, [half_sat_key])

    def _calculate_state_variable_at_patch(self, network, patch_id):
        total_bac = network.get_compartment_value(patch_id, BACTERIUM_EXTRACELLULAR_REPLICATING) + \
                    network.get_compartment_value(patch_id, BACTERIUM_EXTRACELLULAR_DORMANT)
        if not total_bac:
            return 0
        return network.get_compartment_value(patch_id, self._cell_type) * \
               (float(total_bac) / (total_bac + self._parameters[self._half_sat_key]))

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
