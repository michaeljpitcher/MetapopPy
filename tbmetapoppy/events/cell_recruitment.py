from tbmetapoppy.tbpulmonaryenvironment import TBPulmonaryEnvironment
from metapoppy.event import PatchTypeEvent
from parameters import RATE, HALF_SAT


class CellRecruitment(PatchTypeEvent):

    STANDARD_RECRUITMENT = '_standard_recruitment'
    ENHANCED_RECRUITMENT = '_enhanced_recruitment'

    MI_TO_MA_CHEMOKINE_WEIGHT = 'macrophage_infected_to_activated_chemokine_weight'

    def __init__(self, patch_type, cell_type):
        self._cell_type = cell_type
        PatchTypeEvent.__init__(self, patch_type, [cell_type], [], [])

    def _define_parameter_keys(self):
        raise NotImplementedError

    def _calculate_state_variable_at_patch(self, network, patch_id):
        raise NotImplementedError

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._cell_type: 1})


class StandardCellRecruitmentLung(CellRecruitment):
    def __init__(self, cell_type):
        CellRecruitment.__init__(self, TBPulmonaryEnvironment.ALVEOLAR_PATCH, cell_type)
        self._dependent_patch_attributes += [TBPulmonaryEnvironment.PERFUSION]

    def _define_parameter_keys(self):
        rp_key = self._cell_type + CellRecruitment.STANDARD_RECRUITMENT + "_" + self._patch_type + RATE
        return rp_key, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_attribute_value(patch_id, TBPulmonaryEnvironment.PERFUSION)


class EnhancedCellRecruitmentLung(CellRecruitment):
    def __init__(self, cell_recruited):
        self._half_sat_key = None
        self._weight_key = None
        CellRecruitment.__init__(self, TBPulmonaryEnvironment.ALVEOLAR_PATCH, cell_recruited)
        self._dependent_compartments += [TBPulmonaryEnvironment.MACROPHAGE_INFECTED,
                                         TBPulmonaryEnvironment.MACROPHAGE_ACTIVATED]

    def _define_parameter_keys(self):
        rp_key = self._cell_type + CellRecruitment.ENHANCED_RECRUITMENT + "_" + self._patch_type + RATE
        self._half_sat_key = self._cell_type + CellRecruitment.ENHANCED_RECRUITMENT + "_" + self._patch_type + HALF_SAT
        self._weight_key = CellRecruitment.MI_TO_MA_CHEMOKINE_WEIGHT
        return rp_key, [self._half_sat_key, self._weight_key]

    def _calculate_state_variable_at_patch(self, network, patch_id):
        # Get compartment value (will return sum of all values if enhancer is a list of compartments)
        mi = network.get_compartment_value(patch_id, TBPulmonaryEnvironment.MACROPHAGE_INFECTED)
        ma = network.get_compartment_value(patch_id, TBPulmonaryEnvironment.MACROPHAGE_ACTIVATED)
        if not (mi + ma):
            return 0
        ma_and_mi = ma + self._parameters[self._weight_key] * mi
        return network.get_attribute_value(patch_id, TBPulmonaryEnvironment.PERFUSION) * \
               (float(ma_and_mi) / (ma_and_mi + self._parameters[self._half_sat_key]))


class StandardCellRecruitmentLymph(CellRecruitment):
    def __init__(self, cell_type):
        CellRecruitment.__init__(self, TBPulmonaryEnvironment.LYMPH_PATCH, cell_type)

    def _define_parameter_keys(self):
        rp_key = self._cell_type + CellRecruitment.STANDARD_RECRUITMENT + "_" + self._patch_type + RATE
        return rp_key, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return 1


class EnhancedCellRecruitmentLymph(CellRecruitment):
    def __init__(self, cell_type):
        self._half_sat_key = None
        self._weight_key = None
        CellRecruitment.__init__(self, TBPulmonaryEnvironment.LYMPH_PATCH, cell_type)
        self._dependent_compartments += [TBPulmonaryEnvironment.MACROPHAGE_INFECTED,
                                         TBPulmonaryEnvironment.MACROPHAGE_ACTIVATED]

    def _define_parameter_keys(self):
        rp_key = self._cell_type + CellRecruitment.ENHANCED_RECRUITMENT + "_" + self._patch_type + RATE
        self._half_sat_key = self._cell_type + CellRecruitment.ENHANCED_RECRUITMENT + "_" + self._patch_type + HALF_SAT
        self._weight_key = CellRecruitment.MI_TO_MA_CHEMOKINE_WEIGHT
        return rp_key, [self._half_sat_key, self._weight_key]

    def _calculate_state_variable_at_patch(self, network, patch_id):
        # Get compartment value (will return sum of all values if enhancer is a list of compartments)
        mi = network.get_compartment_value(patch_id, TBPulmonaryEnvironment.MACROPHAGE_INFECTED)
        ma = network.get_compartment_value(patch_id, TBPulmonaryEnvironment.MACROPHAGE_ACTIVATED)
        if not (mi + ma):
            return 0
        ma_and_mi = ma + self._parameters[self._weight_key] * mi
        return float(ma_and_mi) / (ma_and_mi + self._parameters[self._half_sat_key])


class EnhancedTCellRecruitmentLymph(CellRecruitment):
    def __init__(self):
        self._half_sat_key = None
        CellRecruitment.__init__(self, TBPulmonaryEnvironment.LYMPH_PATCH, TBPulmonaryEnvironment.T_CELL_NAIVE)
        self._dependent_compartments += [TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE]

    def _define_parameter_keys(self):
        rp_key = self._cell_type + CellRecruitment.ENHANCED_RECRUITMENT + "_" + self._patch_type + RATE
        self._half_sat_key = self._cell_type + CellRecruitment.ENHANCED_RECRUITMENT + "_" + self._patch_type + HALF_SAT
        return rp_key, [self._half_sat_key]

    def _calculate_state_variable_at_patch(self, network, patch_id):
        # Get compartment value (will return sum of all values if enhancer is a list of compartments)
        dcm = network.get_compartment_value(patch_id, TBPulmonaryEnvironment.DENDRITIC_CELL_MATURE)
        if not dcm:
            return 0
        return float(dcm) / (dcm + self._parameters[self._half_sat_key])
