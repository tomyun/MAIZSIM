from .util import VaporPressure
from ..timer import Timer

import copy

#TODO to be predicted using altitude
ATM = 100 # kPa

class Weather:
    def __init__(self, time, PFD, T_air, CO2, RH, wind, P_air=ATM, day_length=None):
        self.setup(time, PFD, T_air, CO2, RH, wind, P_air, day_length)

    def setup(self, time, PFD, T_air, CO2, RH, wind, P_air=ATM, day_length=None):
        self.time = time
        self.PFD = PFD
        self.CO2 = CO2 # ppm
        self.RH = RH # 0~1
        self.T_air = T_air # C
        self.wind = wind # meters s-1
        self.P_air = P_air # kPa
        self.day_length = day_length

    @classmethod
    def from_2DSOIL(cls, T, W):
        #w.jday = W.jday.item()
        #w.time = T.time - W.jday
        self.time = Timer.datetime_from_julian_day(T.time)
        i = T.itime - 1
        self.setup(
            time=time,
            PFD=W.par[i]*4.6, # conversion from PAR in Wm-2 to umol s-1 m-2
            #solRad=W.wattsm[i].item(), # conversion from Wm-2 to J m-2 in one hour Total Radiation incident at soil surface
            T_air=W.tair[i].item(),
            CO2=W.co2.item(),
            RH=np.clip(VaporPressure.relative_humidity(W.vpd[i]), 0.1, 1.0), # % to 0~1
            wind=W.wind*(1000/3600), # conversion from km hr-1 to m s-1
            #rain=W.rint[i].item(), #TODO take account precipitation
            day_length=W.daylng.item(),
        )

    def copy(self):
        return copy.copy(self)

    @property
    def VPD(self):
        return VaporPressure.deficit(self.T_air, self.RH)
