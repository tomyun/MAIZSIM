import numpy as np

class VaporPressure:
    # Campbell and Norman (1998), p 41 Saturation vapor pressure in kPa
    a = 0.611 # kPa
    b = 17.502 # C
    c = 240.97 # C

    #FIXME August-Roche-Magnus formula gives slightly different parameters
    # https://en.wikipedia.org/wiki/Clausius–Clapeyron_relation
    #a = 0.61094 # kPa
    #b = 17.625 # C
    #c = 243.04 # C

    @classmethod
    def saturation(cls, T):
        a, b, c = cls.a, cls.b, cls.c
        return a*np.exp((b*T)/(c+T))

    @classmethod
    def ambient(cls, T, RH):
        es = cls.saturation(T)
        return es * RH

    @classmethod
    def deficit(cls, T, RH):
        es = cls.saturation(T)
        return es * (1 - RH)

    @classmethod
    def relative_humidity(cls, T, VPD):
        es = cls.saturation(T)
        return 1 - VPD / es

    # slope of the sat vapor pressure curve: first order derivative of Es with respect to T
    @classmethod
    def curve_slope(cls, T, P):
        es = cls.saturation(T)
        b, c = cls.b, cls.c
        slope = es * (b*c)/(c+T)**2 / P
        return slope
