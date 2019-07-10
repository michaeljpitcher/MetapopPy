from metapoppy.event import *
from ..tbpulmonaryenvironment import TBPulmonaryEnvironment


class CaseumLiquefy(Event):

    RATE_KEY = 'caseum_liquify_rate'

    def __init__(self):
        Event.__init__(self, [TBPulmonaryEnvironment.MACROPHAGE_ACTIVATED, TBPulmonaryEnvironment.SOLID_CASEUM], [], [])

    def _define_parameter_keys(self):
        return CaseumLiquefy.RATE_KEY, []

    def _calculate_state_variable_at_patch(self, network, patch_id):
        return network.get_compartment_value(patch_id, TBPulmonaryEnvironment.SOLID_CASEUM) * \
               network.get_compartment_value(patch_id, TBPulmonaryEnvironment.MACROPHAGE_ACTIVATED)

    def perform(self, network, patch_id):
        network.update_patch(patch_id, {TBPulmonaryEnvironment.SOLID_CASEUM:-1,
                                        TBPulmonaryEnvironment.LIQUEFIED_CASEUM:1})
