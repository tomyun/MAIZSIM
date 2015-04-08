import numpy as np

class Tracker(object):
    def __init__(self, *args):
        self.reset()
        self.setup(*args)

    def reset(self):
        self._values = []

    def setup(self, *args):
        pass

    def calc(self, T):
        return T

    def update(self, T, dt=1.):
        self._values.append(self.calc(T) * dt)

    @property
    def rate(self):
        return np.mean(self._values)


class ThermalFunc(Tracker):
    @property
    def rate(self):
        return np.sum(self._values)


class BetaFunc(ThermalFunc):
    def setup(self, R_max, T_opt=32.1, T_max=43.7, T_min=0.):
        self.R_max = R_max
        self.T_opt = T_opt * 1.
        self.T_max = T_max * 1.
        self.T_min = T_min * 1.

        # beta function, See Yin et al. (1995), Ag For Meteorol., Yan and Hunt (1999) AnnBot, SK
        self.beta = 1.0
        self.alpha = self.beta * (T_opt - T_min) / (T_max - T_opt)

    def calc(self, T):
        T_opt = self.T_opt
        T_max = self.T_max
        T_min = self.T_min

        if not T_min < T < T_max:
            return 0.
        if not T_min < T_opt < T_maxx:
            return 0.

        f = (T - T_min) / (T_opt - T_min)
        g = (T_max - T) / (T_max - T_opt)
        return self.R_max * f**self.alpha * g**self.beta


class GrowingDegreeDays(ThermalFunc):
    # GDD model with base 8. See Birch et al. (2003) Eu J Agron
    #T_opt = 30.0
    def setup(self, T_base=8.0, T_opt=34.0):
        self.T_base = T_base * 1.
        self.T_opt = T_opt * 1.

    def calc(self, T):
        return min(T, self.T_opt) - self.T_base


#TODO merge vegetative/reproductive versions
class GeneralThermalIndex(ThermalFunc):
    pass


class VegetativeGeneralThermalIndex(ThermalFunc):
    def setup(self, b=0.043177, c=-0.000894):
        self.b = b
        self.c = c

    def calc(self, T):
        return self.b*T**2 + self.c*T**3


class ReproductiveGeneralThermalIndex(ThermalFunc):
    # General Thermal Index, Stewart et al. (1998)
    # Phenological temperature response of maize. Agron. J. 90: 73-79.
    def setup(self, a=5.3581, b=0.011178):
        self.a = a
        self.b = b

    def calc(self, T):
        #T_opt = 32.2
        #a = 1 - 0.6667 * T/T_opt
        return self.a + self.b*T**2

# note it's Tracker, not ThermalFunc
class LeafInductionRate(Tracker):
    def setup(self, temperature_tracker, juvenile_leaves, day_length=None):
        self.temperature_tracker = tempearature_tracker
        self.juvenile_leaves = juvenile_leaves
        self.day_length = day_length

    def calc(self, T):
        #HACK use mean temperature tracker for induction period
        temperature_tracker.update(T)
        T = temperature_tracker.rate

        # effect of photoperiod and temperature on leaf no. used as Grant (1989)
        # Added back the temperature effect on leaf number and revised the algorithm to accumulate addLeafNo to totLeafNo.
        # Changed to respond to mean growing season temperature up to this point.
        # This has little mechanistic basis. Needs improvements. SK 1-19-12

        by_temperature = max(0., 13.6 - 1.89*T + 0.081*T**2 - 0.001*T**3)

        if self.day_length is None:
            by_photo_period = 0.
        else:
            by_photo_period = max(0., 0.1 * (self.juvenile_leaves - 10) * (self.day_length - 12.5))

        return by_temperature + by_photo_period
