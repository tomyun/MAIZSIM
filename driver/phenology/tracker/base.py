import numpy as np

class Tracker(object):
    def __init__(self, timestep=1, *args):
        self.timestep = timestep
        self.reset()
        self.setup(*args)

    def reset(self):
        self._values = []

    def setup(self, *args):
        pass

    def calc(self, T):
        return T

    def set_initial_value(self, value):
        self._values = [value]

    def update(self, T):
        self._values.append(self.calc(T) * self.timestep)

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
