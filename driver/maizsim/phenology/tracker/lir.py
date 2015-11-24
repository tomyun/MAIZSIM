from .base import Tracker

import numpy as np

# note it's Tracker, not Accumulator
class LeafInductionRate(Tracker):
    def setup(self, pheno, temperature_recorder, juvenile_leaves):
        self.pheno = pheno
        self.temperature_recorder = temperature_recorder
        self.juvenile_leaves = juvenile_leaves

    def calc(self, T):
        #TODO implement on_first_update() interface?
        #HACK use mean temperature tracker for induction period
        #TODO use chained methods?
        self.temperature_recorder.update(T)
        T = self.temperature_recorder.rate

        # effect of photoperiod and temperature on leaf no. used as Grant (1989)
        # Added back the temperature effect on leaf number and revised the algorithm to accumulate addLeafNo to totLeafNo.
        # Changed to respond to mean growing season temperature up to this point.
        # This has little mechanistic basis. Needs improvements. SK 1-19-12

        by_temperature = max(0., 13.6 - 1.89*T + 0.081*T**2 - 0.001*T**3)

        day_length = self.pheno.plant.weather.day_length
        if day_length is None:
            by_photo_period = 0.
        else:
            by_photo_period = max(0., 0.1 * (self.juvenile_leaves - 10) * (day_length - 12.5))

        return by_temperature + by_photo_period

    @property
    def rate(self):
        #HACK prevent warnings on nan due to empty _values
        return np.nan if self.empty() else super().rate
