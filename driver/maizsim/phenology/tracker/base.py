import numpy as np

#HACK assumes fixed timestep (probably hourly)
#TODO support dynamic timestep? ever needed?
class Tracker:
    def __init__(self, **kwargs):
        self.reset()
        self.setup(**kwargs)

    def reset(self):
        self._values = []
        return self

    def setup(self, **kwargs):
        pass

    def use_timestep(self, timestep):
        self.timestep = timestep
        return self

    def calc(self, T):
        return T

    def set_initial_value(self, value):
        self._values = [value]
        return self

    def update(self, T):
        self._values.append(self.calc(T) * self.timestep)
        return self

    @property
    def count(self):
        return len(self._values)

    def empty(self):
        return self.count == 0

    @property
    def rate(self):
        return np.mean(self._values)

    @property
    def period(self):
        return self.count * self.timestep


class TemperatureTracker(Tracker):
    def use_timestep(self, timestep):
        #HACK temperature trackers always have timestep of 1
        return super().use_timestep(1)


class Accumulator(Tracker):
    @property
    def rate(self):
        return np.sum(self._values)
