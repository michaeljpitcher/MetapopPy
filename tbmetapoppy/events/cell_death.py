from metapoppy.event import Event
from ..tbcompartments import *
import numpy


class CellDeath(Event):
    def __init__(self, death_rate_key, dying_compartment, additional_parameter_keys=None):
        self._dying_compartment = dying_compartment
        if dying_compartment in INTERNAL_BACTERIA_FOR_CELL:
            self._internal_compartment = INTERNAL_BACTERIA_FOR_CELL[dying_compartment]
        else:
            self._internal_compartment = None
        Event.__init__(self, death_rate_key, additional_parameter_keys)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._dying_compartment)

    def perform(self, network, patch_id):
        changes = {self._dying_compartment: -1}
        if self._internal_compartment:
            bac_to_release = int(round(float(
                network.get_compartment_value(patch_id, self._internal_compartment)) /
                                       network.get_compartment_value(patch_id, self._dying_compartment)))
            changes[self._internal_compartment] = bac_to_release * -1
            changes[BACTERIUM_EXTRACELLULAR_DORMANT] = bac_to_release
        network.update_patch(patch_id, changes)


class MacrophageBursting(CellDeath):
    def __init__(self, death_rate_key, sigmoid_key, capacity_key):
        self._sigmoid_key = sigmoid_key
        self._capacity_key = capacity_key
        CellDeath.__init__(self, death_rate_key, MACROPHAGE_INFECTED, [sigmoid_key, capacity_key])

    def _calculate_state_variable_at_patch(self, network, patch_id):
        bac = network.get_compartment_value(patch_id, BACTERIUM_INTRACELLULAR_MACROPHAGE)
        if not bac:
            return 0
        mac = network.get_compartment_value(patch_id, self._dying_compartment)
        sig = self._parameters[self._sigmoid_key]
        cap = self._parameters[self._capacity_key]
        return mac * ((float(bac) ** sig) / (bac ** sig + ((cap * mac) ** sig)))


class TCellDestroysMacrophage(CellDeath):
    def __init__(self, death_rate_key, half_sat_key):
        self._half_sat_key = half_sat_key
        CellDeath.__init__(self, death_rate_key, MACROPHAGE_INFECTED, [half_sat_key])

    def _calculate_state_variable_at_patch(self, network, patch_id):
        tcell = network.get_compartment_value(patch_id, T_CELL_ACTIVATED)
        if not tcell:
            return 0
        mac = network.get_compartment_value(patch_id, self._dying_compartment)
        return mac * (float(tcell) / (tcell + self._parameters[self._half_sat_key]))


class MacrophageDestroysBacterium(Event):
    def __init__(self, death_rate_key, macrophage_type, half_sat_key):
        self._macrophage_type = macrophage_type
        self._half_sat_key = half_sat_key
        Event.__init__(self, death_rate_key, [half_sat_key])

    def _calculate_state_variable_at_patch(self, network, patch_id):
        mac = network.get_compartment_value(patch_id, self._macrophage_type)
        bac = network.get_compartment_value(patch_id, BACTERIUM_EXTRACELLULAR_REPLICATING) + \
               network.get_compartment_value(patch_id, BACTERIUM_EXTRACELLULAR_DORMANT)
        if not mac or not bac:
            return 0
        return mac * float(bac) / (bac + self._parameters[self._half_sat_key])

    def perform(self, network, patch_id):
        replicating = network.get_compartment_value(patch_id, BACTERIUM_EXTRACELLULAR_REPLICATING)
        dormant = network.get_compartment_value(patch_id, BACTERIUM_EXTRACELLULAR_DORMANT)
        total_bacteria = replicating + dormant

        prob = numpy.array([replicating, dormant], dtype=numpy.dtype('float')) / total_bacteria
        bacteria_type_chosen = numpy.random.choice([BACTERIUM_EXTRACELLULAR_REPLICATING,
                                                    BACTERIUM_EXTRACELLULAR_DORMANT],
                                                   p=prob)

        network.update_patch(patch_id, {bacteria_type_chosen: -1})