from .base import Accumulator

#TODO merge vegetative/reproductive versions
class GeneralThermalIndex(Accumulator):
    pass


class VegetativeGeneralThermalIndex(Accumulator):
    def setup(self, b=0.043177, c=-0.000894):
        self.b = b
        self.c = c

    def calc(self, T):
        return self.b*T**2 + self.c*T**3


class ReproductiveGeneralThermalIndex(Accumulator):
    # General Thermal Index, Stewart et al. (1998)
    # Phenological temperature response of maize. Agron. J. 90: 73-79.
    def setup(self, a=5.3581, b=0.011178):
        self.a = a
        self.b = b

    def calc(self, T):
        #T_opt = 32.2
        #a = 1 - 0.6667 * T/T_opt
        return self.a + self.b*T**2
