from tbmodel import *


class TBDynamicsWithImmuneDrop(TBDynamics):
    # TODO - this does NOT WORK

    RECRUITMENT_DROP_PERCENTAGE = 'recruitment_drop_percentage'
    RECRUITMENT_DROP_INTERVAL = 'recruitment_drop_interval'

    def __init__(self, network_config):
        TBDynamics.__init__(self, network_config)

    def setUp(self, params):
        TBDynamics.setUp(self, params)

        # Post event rate drops
        drop_percent = params[TBDynamicsWithImmuneDrop.RECRUITMENT_DROP_PERCENTAGE]

        current_mr_lung_rate = self._parameters[self._lung_recruit_keys[TBPulmonaryEnvironment.MACROPHAGE_RESTING]]
        current_mr_lymph_rate = self._parameters[self._lymph_recruit_keys[TBPulmonaryEnvironment.MACROPHAGE_RESTING]]
        current_di_rate = self._parameters[self._lung_recruit_keys[TBPulmonaryEnvironment.DENDRITIC_CELL_IMMATURE]]
        current_tn_rate = self._parameters[self._lymph_recruit_keys[TBPulmonaryEnvironment.T_CELL_NAIVE]]

        def drop_rec_rate(rates):
            mr_lung_rate, mr_lymph_rate, di_rate, tn_rate = rates

            self.update_parameter(self._lung_recruit_keys[TBPulmonaryEnvironment.MACROPHAGE_RESTING],
                                  mr_lung_rate)
            self.update_parameter(self._lymph_recruit_keys[TBPulmonaryEnvironment.MACROPHAGE_RESTING],
                                  mr_lymph_rate)
            self.update_parameter(self._lung_recruit_keys[TBPulmonaryEnvironment.DENDRITIC_CELL_IMMATURE],
                                  di_rate)
            self.update_parameter(self._lymph_recruit_keys[TBPulmonaryEnvironment.T_CELL_NAIVE],
                                  tn_rate)

        drop_interval = params[TBDynamicsWithImmuneDrop.RECRUITMENT_DROP_INTERVAL]
        times = [self._start_time + (n * drop_interval) for n in range(1, int(self._max_time/drop_interval)+1)]
        for t in times:
            current_mr_lung_rate = current_mr_lung_rate * (1-drop_percent)
            current_mr_lymph_rate = current_mr_lymph_rate * (1-drop_percent)
            current_di_rate = current_di_rate * (1-drop_percent)
            current_tn_rate = current_tn_rate * (1-drop_percent)
            self.post_event(t, lambda a: drop_rec_rate(a),
                            [current_mr_lung_rate, current_mr_lymph_rate, current_di_rate, current_tn_rate])
