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

    def setup(self, *args):
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

    def empty(self):
        return len(self._values) == 0

    @property
    def rate(self):
        return np.mean(self._values)

    @property
    def period(self):
        return len(self._values) * self.timestep


class Accumulator(Tracker):
    @property
    def rate(self):
        return np.sum(self._values)
