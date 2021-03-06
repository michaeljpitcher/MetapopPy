from ..tbpulmonaryenvironment import *
from metapoppy.event import Event
import numpy
from parameters import HALF_SAT, RATE


class CellIngestBacterium(Event):

    INGEST_BACTERIUM = '_ingest_bacterium'
    INFECTION_PROBABILITY = '_infection_probability'

    def __init__(self, cell_type):
        self._cell_type = cell_type
        if self._cell_type in TBPulmonaryEnvironment.INFECTED_CELL:
            self._infected_cell_type = TBPulmonaryEnvironment.INFECTED_CELL[self._cell_type]
            self._internalised_bac_type = TBPulmonaryEnvironment.INTERNAL_BACTERIA_FOR_CELL[self._infected_cell_type]
        else:
            self._infected_cell_type = self._internalised_bac_type = None
        self._half_sat_key = None
        Event.__init__(self, [self._cell_type, TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT,
                              TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_REPLICATING], [], [])

    def _define_parameter_keys(self):
        rp_key = self._cell_type + CellIngestBacterium.INGEST_BACTERIUM + RATE
        self._half_sat_key = self._cell_type + CellIngestBacterium.INGEST_BACTERIUM + HALF_SAT
        self._infection_prob_key = self._cell_type + CellIngestBacterium.INFECTION_PROBABILITY
        return rp_key, [self._half_sat_key, self._infection_prob_key]

    def _calculate_state_variable_at_patch(self, network, patch_id):
        total_bac = network.get_compartment_value(patch_id, [TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_REPLICATING,
                                                             TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT])
        if not total_bac:
            return 0
        return network.get_compartment_value(patch_id, self._cell_type) * \
           (float(total_bac) / (total_bac + self._parameters[self._half_sat_key]))

    def perform(self, network, patch_id):
        replicating = network.get_compartment_value(patch_id,
                                                    TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_REPLICATING)
        dormant = network.get_compartment_value(patch_id, TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT)
        total_bacteria = replicating + dormant

        r = numpy.random.random() * total_bacteria
        if r < replicating:
            bacteria_type_chosen = TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_REPLICATING
        else:
            bacteria_type_chosen = TBPulmonaryEnvironment.BACTERIUM_EXTRACELLULAR_DORMANT

        # prob = numpy.array([replicating, dormant], dtype=numpy.dtype('float')) / total_bacteria
        # bacteria_type_chosen = numpy.random.choice([TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING,
        #                                             TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT],
        #                                        p=prob)

        r2 = numpy.random.random()
        if r2 < self._parameters[self._infection_prob_key]:
            # Infected
            changes = {self._cell_type: -1, self._infected_cell_type: 1, self._internalised_bac_type: 1,
                       bacteria_type_chosen: -1}
        else:
            # Bacterium destroyed
            changes = {bacteria_type_chosen: -1}
        network.update_patch(patch_id, changes)
