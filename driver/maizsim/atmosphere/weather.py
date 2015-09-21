from .util import VaporPressure
from ..timer import Timer

import numpy as np
import copy

#TODO to be predicted using altitude
ATM = 100 # kPa

class Weather:
    def __init__(self):
        self.reset()

    def reset(self):
        self.time = None
        self.PFD = None # umol m-2 s-1
        self.sol_rad = None # J m-2
        self.CO2 = None # ppm
        self.RH = None # 0~1
        self.T_air = None # C
        self.wind = None # meters s-1
        self.P_air = None # kPa
        self.day_length = None

    @classmethod
    def from_2DSOIL(cls, T, W):
        return Weather().update_from_2DSOIL(T, W)

    def update_from_2DSOIL(self, T, W):
        #w.jday = W.jday.item()
        #w.time = T.time - W.jday
        i = T.itime - 1
        self.time = Timer.datetime_from_julian_day(T.time)
        self.PFD = W.par[i]*4.6 # conversion from PAR in Wm-2 to umol s-1 m-2
        self.sol_rad = W.wattsm[i].item() # conversion from Wm-2 to J m-2 in one hour Total Radiation incident at soil surface
        self.T_air = W.tair[i].item()
        self.CO2 = W.co2.item()
        self.RH = np.clip(VaporPressure.relative_humidity(self.T_air, W.vpd[i]), 0.1, 1.0) # % to 0~1
        self.wind = W.wind*(1000/3600) # conversion from km hr-1 to m s-1
        #self.rain = W.rint[i].item() #TODO take account precipitation
        #HACK use predefined value
        self.P_air = ATM # kPa
        self.day_length = W.daylng.item()
        return self

    def copy(self):
        return copy.copy(self)

    @property
    def VPD(self):
        return VaporPressure.deficit(self.T_air, self.RH)
