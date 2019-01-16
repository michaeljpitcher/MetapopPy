from metapoppy.event import Event
from ..tbpulmonarynetwork import *
import numpy
from parameters import RATE, HALF_SAT, INTRACELLULAR_REPLICATION_SIGMOID, MACROPHAGE_CAPACITY


class CellDeath(Event):
    DEATH = '_death'

    def __init__(self, dying_compartment):
        self._dying_compartment = dying_compartment
        assert dying_compartment not in TBPulmonaryNetwork.INTERNAL_BACTERIA_FOR_CELL, \
            "Standard cell death not appropriate for {0} as it has internal bacteria".format(dying_compartment)
        Event.__init__(self, [self._dying_compartment], [])

    def _define_parameter_keys(self):
        return self._dying_compartment + CellDeath.DEATH + RATE, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, self._dying_compartment)

    def perform(self, network, patch_id):
        changes = {self._dying_compartment: -1}
        network.update_patch(patch_id, changes)


class InfectedCellDeath(CellDeath):
    PERCENT_BACTERIA_DESTROYED = '_percentage_bacteria_destroyed'

    def __init__(self, dying_compartment):
        self._dying_compartment = dying_compartment
        self._internal_bacteria = TBPulmonaryNetwork.INTERNAL_BACTERIA_FOR_CELL[dying_compartment]
        self._bac_percent_to_destroy_key = None
        Event.__init__(self, [self._dying_compartment], [])

    def _define_parameter_keys(self):
        self._bac_percent_to_destroy_key = self._dying_compartment + CellDeath.DEATH + \
                                           InfectedCellDeath.PERCENT_BACTERIA_DESTROYED
        return self._dying_compartment + CellDeath.DEATH + RATE, [self._bac_percent_to_destroy_key]

    def perform(self, network, patch_id):
        avg_bac_per_cell = int(round(float(network.get_compartment_value(patch_id, self._internal_bacteria)) / \
                           network.get_compartment_value(patch_id, self._dying_compartment)))
        bac_to_destroy = int(round(avg_bac_per_cell * self._parameters[self._bac_percent_to_destroy_key]))
        bac_to_release = avg_bac_per_cell - bac_to_destroy
        changes = {self._dying_compartment: -1, self._internal_bacteria: -1 * (bac_to_destroy + bac_to_release),
                   TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT: bac_to_release}
        if self._dying_compartment == TBPulmonaryNetwork.MACROPHAGE_INFECTED:
            changes[TBPulmonaryNetwork.CASEUM] = 1
        network.update_patch(patch_id, changes)


class MacrophageBursting(InfectedCellDeath):

    BURSTING = 'macrophage_bursting'

    def __init__(self):
        InfectedCellDeath.__init__(self, TBPulmonaryNetwork.MACROPHAGE_INFECTED)
        self.dependent_compartments += [TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE]

    def _define_parameter_keys(self):
        self._bac_percent_to_destroy_key = MacrophageBursting.BURSTING + InfectedCellDeath.PERCENT_BACTERIA_DESTROYED
        return MacrophageBursting.BURSTING + RATE, \
               [self._bac_percent_to_destroy_key, MACROPHAGE_CAPACITY, INTRACELLULAR_REPLICATION_SIGMOID]

    def _calculate_state_variable_at_patch(self, network, patch_id):
        bac = network.get_compartment_value(patch_id, TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE)
        if not bac:
            return 0
        mac = network.get_compartment_value(patch_id, self._dying_compartment)
        sig = self._parameters[INTRACELLULAR_REPLICATION_SIGMOID]
        cap = self._parameters[MACROPHAGE_CAPACITY]
        return mac * ((float(bac) ** sig) / (bac ** sig + ((cap * mac) ** sig)))


class TCellDestroysMacrophage(InfectedCellDeath):

    T_CELL_DESTROYS_MACROPHAGE = 't_cell_destroys_macrophage'

    def __init__(self):
        InfectedCellDeath.__init__(self, TBPulmonaryNetwork.MACROPHAGE_INFECTED)
        self.dependent_compartments += [TBPulmonaryNetwork.BACTERIUM_INTRACELLULAR_MACROPHAGE]

    def _define_parameter_keys(self):
        self._bac_percent_to_destroy_key = TCellDestroysMacrophage.T_CELL_DESTROYS_MACROPHAGE + \
                                           InfectedCellDeath.PERCENT_BACTERIA_DESTROYED
        self._half_sat_key = TCellDestroysMacrophage.T_CELL_DESTROYS_MACROPHAGE + HALF_SAT
        return TCellDestroysMacrophage.T_CELL_DESTROYS_MACROPHAGE + RATE, \
               [self._bac_percent_to_destroy_key, self._half_sat_key]

    def _calculate_state_variable_at_patch(self, network, patch_id):
        tcell = network.get_compartment_value(patch_id, TBPulmonaryNetwork.T_CELL_ACTIVATED)
        if not tcell:
            return 0
        mac = network.get_compartment_value(patch_id, self._dying_compartment)
        return mac * (float(tcell) / (tcell + self._parameters[self._half_sat_key]))


class MacrophageDestroysBacterium(Event):

    DESTROYS_BACTERIUM = '_destroys_bacterium'

    def __init__(self, macrophage_type):
        self._macrophage_type = macrophage_type
        self._half_sat_key = None
        Event.__init__(self, [TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING,
                              TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT,
                              self._macrophage_type], [])

    def _define_parameter_keys(self):
        rp_key = self._macrophage_type + MacrophageDestroysBacterium.DESTROYS_BACTERIUM + RATE
        self._half_sat_key = self._macrophage_type + MacrophageDestroysBacterium.DESTROYS_BACTERIUM + HALF_SAT
        return rp_key, [self._half_sat_key]

    def _calculate_state_variable_at_patch(self, network, patch_id):
        mac = network.get_compartment_value(patch_id, self._macrophage_type)
        bac = network.get_compartment_value(patch_id, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING) + \
               network.get_compartment_value(patch_id, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT)
        if not mac or not bac:
            return 0
        return mac * float(bac) / (bac + self._parameters[self._half_sat_key])

    def perform(self, network, patch_id):
        replicating = network.get_compartment_value(patch_id, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING)
        dormant = network.get_compartment_value(patch_id, TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT)
        total_bacteria = replicating + dormant

        prob = numpy.array([replicating, dormant], dtype=numpy.dtype('float')) / total_bacteria
        bacteria_type_chosen = numpy.random.choice([TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_REPLICATING,
                                                    TBPulmonaryNetwork.BACTERIUM_EXTRACELLULAR_DORMANT],
                                                   p=prob)

        network.update_patch(patch_id, {bacteria_type_chosen: -1})