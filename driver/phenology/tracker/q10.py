from .base import Accumulator

class Q10Func(Accumulator):
    def setup(self, T_opt, Q10=2.0):
        self.T_opt = T_opt
        self.Q10 = Q10

    def calc(self, T):
        return self.Q10 ** ((T - self.T_opt) / 10)
