from .base import Accumulator

class GrowingDegreeDays(Accumulator):
    # GDD model with base 8. See Birch et al. (2003) Eu J Agron
    #T_opt = 30.0
    def setup(self, T_base=8.0, T_opt=34.0, T_max=None):
        self.T_base = T_base
        self.T_opt = T_opt
        self.T_max = T_max

    def calc(self, T):
        if self.T_opt is not None:
            T = min(T, self.T_opt)
        if self.T_max is not None:
            T = self.T_base if T >= self.T_max else T
        return T - self.T_base
