from .base import Accumulator

class BetaFunc(Accumulator):
    def setup(self, R_max, T_opt=32.1, T_max=43.7, T_min=0.):
        self.R_max = R_max
        self.T_opt = T_opt
        self.T_max = T_max
        self.T_min = T_min

        # beta function, See Yin et al. (1995), Ag For Meteorol., Yan and Hunt (1999) AnnBot, SK
        self.beta = 1.0
        self.alpha = self.beta * (T_opt - T_min) / (T_max - T_opt)

    def calc(self, T):
        T_opt = self.T_opt
        T_max = self.T_max
        T_min = self.T_min

        if not T_min < T < T_max:
            return 0
        if not T_min < T_opt < T_max:
            return 0

        f = (T - T_min) / (T_opt - T_min)
        g = (T_max - T) / (T_max - T_opt)
        return self.R_max * f**self.alpha * g**self.beta
