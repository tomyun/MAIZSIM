from .util import VaporPressure

import copy

class Weather:
    def __init__(self, PFD, T_air, CO2, RH, wind, P_air):
        self.setup(PFD, T_air, CO2, RH, wind, P_air)

    def setup(self, PFD, T_air, CO2, RH, wind, P_air):
        self.PFD = PFD
        self.CO2 = CO2 # ppm
        self.RH = np.clip(RH, 10, 100) / 100. # %
        self.T_air = T_air # C
        self.wind = wind # meters s-1
        self.P_air = P_air # kPa

    def copy(self):
        return copy.copy(self)

    @property
    def VPD(self):
        return VaporPressure.deficit(self.T_air, self.RH)
