from tbmodel import *


class TBDynamicsWithHIV(TBDynamics):

    HIV_INITIAL_TIME = 'hiv_initial_time'
    HIV_DROP = 'hiv_recruit_drop'

    def __init__(self, network_config):
        TBDynamics.__init__(self, network_config)

    def setUp(self, params):
        TBDynamics.setUp(self, params)

        assert 0 <= params[TBDynamicsWithHIV.HIV_DROP] <= 1.0, "HIV drop must be % (0-1)"

        # Get params
        hiv_initial_time = params[TBDynamicsWithHIV.HIV_INITIAL_TIME]
        hiv_drop = params[TBDynamicsWithHIV.HIV_DROP]

        # Calculate new rate
        initial_tn_rate = self._parameters[self._lymph_recruit_keys[TBPulmonaryEnvironment.T_CELL_NAIVE]]
        new_tn_rate = initial_tn_rate * (1-hiv_drop)

        def update_rate(rates):
            self.update_parameter(self._lymph_recruit_keys[TBPulmonaryEnvironment.T_CELL_NAIVE], rates[0])

        self.post_event(self._start_time, lambda a: update_rate(a), [initial_tn_rate])
        self.post_event(hiv_initial_time, lambda a: update_rate(a), [new_tn_rate])
