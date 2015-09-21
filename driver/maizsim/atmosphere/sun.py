# Unit to calculate solar geometry including solar elevation, declination,
#  azimuth etc using TSolar class. Data are hidden. 03/15/00 SK
# - 1st Revision 10/10/00: Changed to dealing only upto the top of the canopy. Radiation transter within the canopy is now a separate module.
# - added functins to calculate global radiation, atmospheric emissivity, etc as in Spitters et al. (1986), 3/18/01, SK
# 24Dec03, SK
# - added a function to calculate day length based on Campbell and Norman (1998) p 170,
# - added overloads to SetVal
# 2Aug04, SK
# - Translated to C++ from Delphi
# - revised some functions according to "An introduction to solar radiaiton" by Iqbal (1983)
# - (To Do) Add algorithms for instantaneous diffuse and direct radiation predictions from daily global solar radiation for given time
# - (To Do) This can be done by first applying sinusoidal model to the daily data to simulate hourly global solar radiation
# - (To Do) Then the division model of diffuse and direct radiations was applied
# - added direct and diffuse separation functions according to Weiss and Norman (1985), 3/16/05

from ..timer import Timer

from numpy import pi, sin, cos, tan, arcsin, arccos, radians, degrees, log, exp, fmax

# conversion factor from W/m2 to PFD (umol m-2 s-1) for PAR waveband (median 550 nm of 400-700 nm) of solar radiation,
# see Campbell and Norman (1994) p 149
# 4.55 is a conversion factor from W to photons for solar radiation, Goudriaan and van Laar (1994)
# some use 4.6 i.e., Amthor 1994, McCree 1981, Challa 1995.
PHOTON_UMOL_PER_J = 4.6

# solar constant, Iqbal (1983)
#FIXME better to be 1361 or 1362 W/m-2?
SOLAR_CONSTANT = 1370.0

class Sun:
    def __init__(self, time, latitude, longitude, altitude=50, transmissivity=0.5, PAR=None):
        self.setup(time, latitude, longitude, altitude, transmissivity, PAR)

    def setup(self, time, latitude, longitude, altitude, transmissivity, PAR):
        #HACK takes account different Julian day conventions (03-01 vs. 01-01)
        self.time = time
        self.day = int(time.strftime('%j'))
        self.hour = time.hour
        self.latitude = latitude # DO NOT convert to radians for consistency
        self.longitude = longitude # leave it as in degrees, used only once for solar noon calculation
        self.altitude = altitude # m from the sea level
        self.transmissivity = transmissivity # atmospheric transmissivity, Goudriaan and van Laar (1994) p 30
        self.PAR = PAR # observed PAR in PPFD (umol m-2 s-1) unit

    #####################
    # Solar Coordinates #
    #####################

    #HACK always use degrees for consistency and easy tracing
    @property
    def declination_angle(self):
        # Goudriaan 1977
        def goudriaan(d):
            g = 2*pi*(d + 10) / 365
            return -23.45 * cos(g)

        # Resenberg, blad, verma 1982
        def resenberg(d):
            g = 2*pi*(d - 172) / 365
            return 23.5 * cos(g)

        # Iqbal (1983) Pg 10 Eqn 1.3.3, and sundesign.com
        def iqbal(d):
            g = 2*pi*(d + 284) / 365
            return 23.45 * sin(g)

        # Campbell and Norman, p168
        def campbell(d):
            a = radians(356.6 + 0.9856*d)
            b = radians(278.97 + 0.9856*d + 1.9165*sin(a))
            r = arcsin(0.39785*sin(b))
            return degrees(r)

        # Spencer equation, Iqbal (1983) Pg 7 Eqn 1.3.1. Most accurate among all
        def spencer(d):
            # gamma: day angle
            g = 2*pi*(d - 1) / 365
            r = 0.006918 - 0.399912*cos(g) + 0.070257*sin(g) - 0.006758*cos(2*g) + 0.000907*sin(2*g) -0.002697*cos(3*g) + 0.00148*sin(3*g)
            return degrees(r)
        #FIXME pascal version of LightEnv uses iqbal()
        return spencer(self.day)


    # LC is longitude correction for Light noon, Wohlfart et al, 2000; Campbell & Norman 1998
    @property
    def longitude_correction(self):
        # standard meridian for pacific time zone is 120 W, Eastern Time zone : 75W
        # LC is positive if local meridian is east of standard meridian, i.e., 76E is east of 75E
        standard_meridian = 120
        #FIXME shouldn't it be negative for West?
        #HACK assume we're on the west
        return -(self.longitude - standard_meridian) / (360 / 24)

    @property
    def solar_noon(self):
        LC = self.longitude_correction
        # epsilon?: convert degrees to radians
        f = radians(279.575 + 0.9856*self.day)
        # calculating Equation of Time
        ET = (-104.7*sin(f) + 596.2*sin(2*f) + 4.3*sin(3*f) - 12.7*sin(4*f) \
              -429.3*cos(f) - 2.0*cos(2*f) + 19.3*cos(3*f)) / (60 * 60)
        return 12 - LC - ET

    def _cos_hour_angle(self, zenith_angle):
        # this value should never become negative because -90 <= latitude <= 90 and -23.45 < decl < 23.45
        #HACK is this really needed for crop models?
        # preventing division by zero for N and S poles
        #denom = fmax(denom, 0.0001)
        # sunrise/sunset hour angle
        #TODO need to deal with lat_bound to prevent tan(90)?
        #lat_bound = radians(68)? radians(85)?
        # cos(h0) at cos(theta_s) = 0 (solar zenith angle = 90 deg == elevation angle = 0 deg)
        #return -tan(self.latitude) * tan(self.declination_angle)
        w_s = radians(zenith_angle)
        p = radians(self.latitude)
        d = radians(self.declination_angle)
        return (cos(w_s) - sin(p) * sin(d)) / (cos(p) * cos(d))

    @property
    def half_day_length(self):
        # from Iqbal (1983) p 16
        def hour_angle_at_horizon():
            c = self._cos_hour_angle(90)
            # in the polar region during the winter, sun does not rise
            if c > 1:
                return 0
            # white nights during the summer in the polar region
            elif c < -1:
                return 180
            else:
                return degrees(arccos(c))
        return hour_angle_at_horizon() / (360 / 24)

    @property
    def day_length(self):
        return self.half_day_length * 2

    @property
    def sunrise(self):
        return self.solar_noon - self.half_day_length

    @property
    def sunset(self):
        return self.solar_noon + self.half_day_length

    @property
    def hour_angle(self):
        return (self.hour - self.solar_noon) * (360 / 24)

    @property
    def elevation_angle(self):
        #FIXME When time gets the same as solarnoon, this function fails. 3/11/01 ??
        h = radians(self.hour_angle)
        p = radians(self.latitude)
        d = radians(self.declination_angle)
        return arcsin(cos(h) * cos(d) * cos(p) + sin(d) * sin(p))

    @property
    def zenith_angle(self):
        #FIXME need abs()?
        return abs(90 - self.elevation_angle)

    # The solar azimuth angle is the angular distance between due South and the
    # projection of the line of sight to the sun on the ground.
    # View point from south, morning: +, afternoon: -
	# See An introduction to solar radiation by Iqbal (1983) p 15-16
	# Also see http://www.susdesign.com/sunangle/index.html
    @property
    def azimuth_angle(self):
        d = radians(self.declination_angle)
        t_s = radians(self.elevation_angle)
        p = radians(self.latitude)
        r = arccos((sin(d) - sin(t_s) * sin(p)) / (cos(t_s) * cos(p)))
        return degrees(abs(r))

    ###################
    # Solar Radiation #
    ###################

    @property
    def solar_radiation(self):
        # Campbell and Norman's global solar radiation, this approach is used here
        return self.directional_solar_radiation + self.diffusive_solar_radiation

    # atmospheric pressure in kPa
    @property
    def atmospheric_pressure(self):
        try:
            # campbell and Norman (1998), p 41
            return 101.3 * exp(-self.altitude / 8200)
        except:
            return 100

    @property
    def optical_air_mass_number(self):
        t_s = fmax(0, radians(self.elevation_angle))
        #FIXME need to do np.fmax(0.0001, sin(t_s))?
        return self.atmospheric_pressure / (101.3 * sin(t_s))

    @property
    #TODO rename to insolation?
    def _solar_radiation(self):
        t_s = fmax(0, radians(self.elevation_angle))
        g = 2*pi*(self.day - 10) / 365
        return SOLAR_CONSTANT * sin(t_s) * (1 + 0.033*cos(g))

    @property
    def directional_solar_radiation(self):
        return self.directional_coeff * self._solar_radiation

    @property
    def diffusive_solar_radiation(self):
        return self.diffusive_coeff * self._solar_radiation

    @property
    def directional_coeff(self):
        # Goudriaan and van Laar's global solar radiation
        def goudriaan(tau):
            #FIXME should be goudriaan() version
            return tau * (1 - self.diffusive_coeff)

        # Takakura (1993), p 5.11
        def takakura(tau, m):
            return tau**m

        # Campbell and Norman (1998), p 173
        def campbell(tau, m):
            return tau**m
        return campbell(self.transmissivity, self.optical_air_mass_number)

    # Fdif: Fraction of diffused light
    @property
    def diffusive_coeff(self):
        # Goudriaan and van Laar's global solar radiation
        def goudriaan(tau):
            # clear sky : 20% diffuse
            if tau >= 0.7:
                return 0.2
            # cloudy sky: 100% diffuse
            elif tau <= 0.3:
                return 1
            # inbetween
            else:
                return 1.6 - 2*tau

        # Takakura (1993), p 5.11
        def takakura(tau, m):
            return (1 - tau**m) / (1 - 1.4*log(tau)) / 2

        # Campbell and Norman (1998), p 173
        def campbell(tau, m):
            return (1 - tau**m) * 0.3
        return campbell(self.transmissivity, self.optical_air_mass_number)

    @property
    def directional_fraction(self):
        return 1 / (1 + self.diffusive_coeff / self.directional_coeff)

    @property
    def diffusive_fraction(self):
        return 1 / (1 + self.directional_coeff / self.diffusive_coeff)

    # PARfr
    @property
    #TODO better naming: extinction? transmitted_fraction?
    def photosynthetic_coeff(self):
        #if self.elevation_angle <= 0:
        #    return 0

        #TODO: implement Weiss and Norman (1985), 3/16/05
        def weiss():
            pass

        # Goudriaan and van Laar (1994)
        def goudriaan(tau):
            # clear sky (tau >= 0.7): 45% is PAR
            if tau >= 0.7:
                return 0.45
            # cloudy sky (<= 0.3): 55% is PAR
            elif tau <= 0.3:
                return 0.55
            else:
                return 0.625 - 0.25*tau
        return goudriaan(self.transmissivity)

    # PARtot: total PAR (umol m-2 s-1) on horizontal surface
    @property
    def photosynthetic_radiation(self):
        if self.PAR is not None:
            return self.PAR
        else:
            return self.solar_radiation * self.photosynthetic_coeff * PHOTON_UMOL_PER_J

    # PARdir
    @property
    def directional_photosynthetic_radiation(self):
        return self.directional_fraction * self.photosynthetic_radiation

    # PARdif
    @property
    def diffusive_photosynthetic_radiation(self):
        return self.diffusive_fraction * self.photosynthetic_radiation
