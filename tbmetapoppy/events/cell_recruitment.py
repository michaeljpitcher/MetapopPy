from tbmetapoppy.pulmonarynetwork import PulmonaryNetwork
from metapoppy.event import PatchTypeEvent


class CellRecruitment(PatchTypeEvent):
    def __init__(self, patch_type, recruitment_rate_key, cell_type, additional_parameters=None):
        self._cell_type = cell_type
        PatchTypeEvent.__init__(self, patch_type, recruitment_rate_key, additional_parameters)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        raise NotImplementedError

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {self._cell_type: 1})


class StandardCellRecruitmentLung(CellRecruitment):
    def __init__(self, recruitment_rate_key, cell_type):
        CellRecruitment.__init__(self, PulmonaryNetwork.ALVEOLAR_PATCH, recruitment_rate_key, cell_type)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_attribute_value(patch_id, PulmonaryNetwork.PERFUSION)


class EnhancedCellRecruitmentLung(CellRecruitment):
    def __init__(self, recruitment_rate_key, cell_type, enhancer, half_sat_key):
        self._half_sat_key = half_sat_key
        self._enhancers = enhancer
        CellRecruitment.__init__(self, PulmonaryNetwork.ALVEOLAR_PATCH, recruitment_rate_key, cell_type, [half_sat_key])

    def _calculate_state_variable_at_patch(self, network, patch_id):
        # Get compartment value (will return sum of all values if enhancer is a list of compartments)
        enhancer_count = network.get_compartment_value(patch_id, self._enhancers)
        if not enhancer_count:
            return 0
        return network.get_attribute_value(patch_id, PulmonaryNetwork.PERFUSION) * (float(enhancer_count) /
                (enhancer_count + self._parameters[self._half_sat_key]))


class StandardCellRecruitmentLymph(CellRecruitment):
    def __init__(self, recruitment_rate_key, cell_type):
        CellRecruitment.__init__(self, PulmonaryNetwork.LYMPH_PATCH, recruitment_rate_key, cell_type)

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return 1


class EnhancedCellRecruitmentLymph(CellRecruitment):
    def __init__(self, recruitment_rate_key, cell_type, enhancer, half_sat_key):
        self._half_sat_key = half_sat_key
        self._enhancers = enhancer
        CellRecruitment.__init__(self, PulmonaryNetwork.LYMPH_PATCH, recruitment_rate_key, cell_type, [half_sat_key])

    def _calculate_state_variable_at_patch(self, network, patch_id):
        # Get compartment value (will return sum of all values if enhancer is a list of compartments)
        enhancer_count = network.get_compartment_value(patch_id, self._enhancers)
        if not enhancer_count:
            return 0
        return float(enhancer_count) / (enhancer_count + self._parameters[self._half_sat_key])
