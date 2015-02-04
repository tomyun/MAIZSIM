# Coupled model of photosynthesis-stomatal conductance-energy balance for a maize leaf
# this unit simulates Maize leaf gas-exchange characteristics
# including photosynthesis, traspiration, boundary and stomatal conductances,
# and leaf temperature based on von Caemmerer (2000) C4 model, BWB stomatal
# conductance (1987) and Energy balance model as described in Campbell and
# Norman (1998) photosynthetic parameters were calibrated with PI3733 from
# SPAR experiments at Beltsville, MD in 2002 Stomatal conductance parameters
# were not calibrated
# Ver 1.0, S.Kim, 11/2002, Originally written in Pascal
# Translated into C++ by S.Kim June 11, 2003 following D.Timlin's translation of C3 model
# IMPORTANT: This model must not be released until validated and published
# 2006-2007 modified to use leaf water potentials to adjust photosynthesis for water stress Y. Yang
# modified 2009 to adjust photosynthesis for nitroge stress Y. Yang

import numpy as np
import scipy.optimize

#FIXME are they parameters?
EPS = 0.97
SBC = 5.6697e-8

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
    def saturation(cls, t):
        a, b, c = cls.a, cls.b, cls.c
        return a*np.exp((b*t)/(c+t))

    @classmethod
    def ambient(cls, t, rh):
        es = cls.saturation(t)
        return es * rh

    @classmethod
    def deficit(cls, t, rh):
        es = cls.saturation(t)
        return es * (1 - rh)

    # slope of the sat vapor pressure curve: first order derivative of Es with respect to T
    @classmethod
    def curve_slope(cls, t, press):
        es = cls.saturation(t)
        b, c = cls.b, cls.c
        slope = es * (b*c)/(c+t)**2 / press
        return slope


class Atmosphere:
    def __init__(self, pfd, tair, co2, rh, wind, press):
        self.setup(pfd, tair, co2, rh, wind, press)

    def setup(self, pfd, tair, co2, rh, wind, press):
        self.pfd = pfd
        par = pfd / 4.55
        # If total solar radiation unavailable, assume NIR the same energy as PAR waveband
        nir = par

        scatt = 0.15
        # times 2 for projected area basis
        self.r_abs = (1 - scatt)*par + 0.15*nir + 2*(EPS * SBC * (tair+273)**4)

        # shortwave radiation (PAR (=0.85) + NIR (=0.15) solar radiation absorptivity of leaves: =~ 0.5
        self.co2 = co2
        self.rh = np.clip(rh, 10, 100) / 100.
        self.tair = tair # C
        self.wind = wind # meters s-1
        self.press = press # kPa


class Stomata:
    def __init__(self, leaf_width):
        self.setup(leaf_width)

    def setup(self, leaf_width):
        # in P. J. Sellers, et al.Science 275, 502 (1997)
        # g0 is b, of which the value for c4 plant is 0.04
        # and g1 is m, of which the value for c4 plant is about 4 YY
        self.g0 = 0.04
        self.g1 = 4.0

        self.gb = 0. # boundary layer conductance
        self.gs = 0. # stomatal conductance

        self.leaf_width = leaf_width / 100. # meters
        self.leafp_effect = 1 # At first assume there is not drought stress, so assign 1 to leafpEffect. Yang 8/20/06

    def update_boundary_layer(self, wind):
        # maize is an amphistomatous species, assume 1:1 (adaxial:abaxial) ratio.
        sr = 1.0
        ratio = (sr + 1)**2 / (sr**2 + 1)

        # characteristic dimension of a leaf, leaf width in m
        d = self.leaf_width * 0.72

        #return 1.42 # total BLC (both sides) for LI6400 leaf chamber
        self.gb = 1.4 * 0.147 * (max(0.1, wind) / d)**0.5 * ratio
        #self.gb = (1.4 * 1.1 * 6.62 * (wind / d)**0.5 * (press / (R * (273.15 + tair)))) # this is an alternative form including a multiplier for conversion from mm s-1 to mol m-2 s-1
        # 1.1 is the factor to convert from heat conductance to water vapor conductance, an avarage between still air and laminar flow (see Table 3.2, HG Jones 2014)
        # 6.62 is for laminar forced convection of air over flat plates on projected area basis
        # when all conversion is done for each surface it becomes close to 0.147 as given in Norman and Campbell
        # multiply by 1.4 for outdoor condition, Campbell and Norman (1998), p109, also see Jones 2014, pg 59 which suggest using 1.5 as this factor.
        # multiply by ratio to get the effective blc (per projected area basis), licor 6400 manual p 1-9
        return self.gb

    # stomatal conductance for water vapor in mol m-2 s-1
    def update_stomata(self, pressure, co2, a_net, rh, tleaf):
        # params
        g0 = self.g0
        g1 = self.g1
        gb = self.gb

        gamma = 10.0
        cs = co2 - (1.37 * a_net / gb) # surface CO2 in mole fraction
        if cs <= gamma:
            cs = gamma + 1

        temp = self._set_leafp_effect(pressure)

        aa = temp * g1 * a_net / cs
        bb = g0 + gb - (temp * g1 * a_net / cs)
        cc = (-rh * gb) - g0
        #hs = QuadSolnUpper(aa, bb, cc)
        hs = max(np.roots([aa, bb, cc]))
        hs = np.clip(hs, 0.3, 1.) # preventing bifurcation

        #es = VaporPressure.saturation(tleaf)
        #ds = (1 - hs) * es # VPD at leaf surface
        ds = VaporPressure.deficit(tleaf, hs)
        tmp = g0 + (g1 *temp * (a_net * hs / cs))
        tmp = max(tmp, g0)

        # this below is an example of how you can write temporary data to a debug window. It can be copied and
        # pasted into excel for plotting. Dennis See above where the CString object is created.
        print("tmp = %f pressure = %f Ds= %f Tleaf = %f Cs = %f Anet = %f hs = %f RH = %f" % (tmp, pressure, ds, tleaf, cs, a_net, hs, rh))
        self.gs = tmp
        return self.gs

    def _set_leafp_effect(self, pressure):
        # pressure - leaf water potential MPa...
        sf = 2.3 # sensitivity parameter Tuzet et al. 2003 Yang
        phyf = -1.2 # reference potential Tuzet et al. 2003 Yang
        self.leafp_effect = (1 + np.exp(sf * phyf)) / (1 + np.exp(sf * (phyf - pressure)))
        return self.leafp_effect

    def total_conductance_h20(self):
        gs = self.gs
        gb = self.gb
        return gs * gb / (gs + gb)

    def boundary_layer_resistance_co2(self):
        return 1.37 / self.gb

    def stomatal_resistance_co2(self):
        return 1.6 / self.gs

    def total_resistance_co2(self):
        return self.boundary_layer_resistance_co2() + self.stomatal_resistance_co2()


class Photosynthesis:
    def __init__(self, stomata, leaf_n_content):
        self.setup()
        self.stomata = stomata
        self.leaf_n_content = leaf_n_content

    def setup(self):
        # activation energy values
        self.eac = 59400.
        self.eao = 36000.

        self.eavp = 75100.
        self.eavc = 55900. # Sage (2002) JXB
        self.eaj = 32800.

        self.hj = 220000.
        self.sj = 702.6

        self.kc25 = 650. # Michaelis constant of rubisco for CO2 of C4 plants (2.5 times that of tobacco), ubar, Von Caemmerer 2000
        self.ko25 = 450. # Michaelis constant of rubisco for O2 (2.5 times C3), mbar
        self.kp25 = 80. # Michaelis constant for PEP caboxylase for CO2

        # Kim et al. (2007), Kim et al. (2006)
        # In von Cammerer (2000), Vpm25=120, Vcm25=60,Jm25=400
        # In Soo et al.(2006), under elevated C5O2, Vpm25=91.9, Vcm25=71.6, Jm25=354.2 YY
        self.vpm25 = 70.
        self.vcm25 = 50.
        self.jm25 = 300.

        # Values in Kim (2006) are for 31C, and the values here are normalized for 25C. SK
        self.rd25 = 2.
        self.ear = 39800.

        #FIXME are they even used?
        #self.beta_aba = 1.48e2 # Tardieu-Davies beta, Dewar (2002) Need the references !?
        #self.delta = -1.0
        #self.alpha_aba = 1.0e-4
        #self.lambda_r = 4.0e-12 # Dewar's email
        #self.lambda_l = 1.0e-12
        #self.k_max = 6.67e-3 # max. xylem conductance (mol m-2 s-1 MPa-1) from root to leaf, Dewar (2002)

    def _light(self, pfd):
        #FIXME make scatt global parameter?
        scatt = 0.15 # leaf reflectance + transmittance
        f = 0.15 #spectral correction

        ia = pfd * (1 - scatt) # absorbed irradiance
        i2 = ia * (1 - f) / 2. # useful light absorbed by PSII
        return i2

    # mesophyll CO2 partial pressure, ubar, one may use the same value as Ci assuming infinite mesohpyle conductance
    def _co2_mesophyll(self, a_net, press, co2, stomata):
        p = press / 100.
        ca = co2 * p # conversion to partial pressure
        cm = ca - a_net * stomata.total_resistance_co2() * p
        return np.clip(cm, 0., 2*ca)

    # Arrhenius equation
    def _temperature_dependence_rate(self, Ea, T, Tb=25.):
        R = 8.314 # universal gas constant (J K-1 mol-1)
        K = 273.
        return np.exp(Ea * (T - Tb) / ((Tb + K) * R * (T + K)))

    def _nitrogen_limited_rate(self, N):
        # in Sinclair and Horie, 1989 Crop sciences, it is 4 and 0.2
        # In J Vos. et al. Field Crop study, 2005, it is 2.9 and 0.25
        # In Lindquist, weed science, 2001, it is 3.689 and 0.5
        s = 2.9 # slope
        N0 = 0.25
        return 2 / (1 + np.exp(-s * (max(N0, N) - N0))) - 1

    def _dark_respiration(self, tleaf):
        return self.rd25 * self._temperature_dependence_rate(self.ear, tleaf)

    def _maximum_electron_transport_rate(self, T, N):
        R = 8.314

        Tb = 25.
        K = 273.
        Tk = T + K
        Tbk = Tb + K

        Sj = self.sj
        Hj = self.hj

        return self.jm25 * self._nitrogen_limited_rate(N) \
                         * self._temperature_dependence_rate(self.eaj, T) \
                         * (1 + np.exp((Sj*Tbk - Hj) / (R*Tbk))) \
                         / (1 + np.exp((Sj*Tk  - Hj) / (R*Tk)))

    # Incident PFD, Air temp in C, CO2 in ppm, RH in percent
    #FIXME remove dependency on leaf parameters (pressure=LWP, tleaf)
    def photosynthesize(self, pfd, press, co2, rh, pressure, tleaf):
        # Light response function parameters
        i2 = self._light(pfd)

        o = 210. # gas units are mbar
        om = o # mesophyll O2 partial pressure

        kp = self.kp25 # T dependence yet to be determined
        kc = self.kc25 * self._temperature_dependence_rate(self.eac, tleaf)
        ko = self.ko25 * self._temperature_dependence_rate(self.eao, tleaf)
        km = kc * (1 + om / ko) # effective M-M constant for Kc in the presence of O2

        rd = self._dark_respiration(tleaf)
        rm = 0.5 * rd

        vpmax = self.vpm25 * self._nitrogen_limited_rate(self.leaf_n_content) * self._temperature_dependence_rate(self.eavp, tleaf)
        vcmax = self.vcm25 * self._nitrogen_limited_rate(self.leaf_n_content) * self._temperature_dependence_rate(self.eavc, tleaf)
        jmax  = self._maximum_electron_transport_rate(tleaf, self.leaf_n_content)

        def c4(a_net, press, co2, rh, tleaf):
            gbs = 0.003 # bundle sheath conductance to CO2, mol m-2 s-1
            #gi = 1.0 # conductance to CO2 from intercelluar to mesophyle, mol m-2 s-1, assumed

            self.stomata.update_stomata(pressure, co2, a_net, rh, tleaf)
            cm = self._co2_mesophyll(a_net, press, co2, self.stomata)

            def enzyme_limited():
                # PEP carboxylation rate, that is the rate of C4 acid generation
                vp1 = (cm * vpmax) / (cm + kp)
                vp2 = vpr = 80. # PEP regeneration limited Vp, value adopted from vC book
                vp = max(min(vp1, vp2), 0)

                #FIXME where should gamma be at?
                # half the reciprocal of rubisco specificity, to account for O2 dependence of CO2 comp point,
                # note that this become the same as that in C3 model when multiplied by [O2]
                #gamma1 = 0.193
                #gamma_star = gamma1 * os
                #gamma = (rd*km + vcmax*gamma_star) / (vcmax - rd)

                # Enzyme limited A (Rubisco or PEP carboxylation
                ac1 = vp + gbs*cm - rm
                #ac1 = max(0, ac1) # prevent Ac1 from being negative Yang 9/26/06
                ac2 = vcmax - rd
                ac = min(ac1, ac2)
                return ac
            ac = enzyme_limited()

            # Light and electron transport limited A mediated by J
            def transport_limited():
                theta = 0.5
                j = min(np.roots([theta, -(i2+jmax), i2*jmax])) # rate of electron transport
                x = 0.4 # Partitioning factor of J, yield maximal J at this value
                aj1 = x*j/2. - rm + gbs*cm
                aj2 = (1 - x)*j/3. - rd
                aj = min(aj1, aj2)
                return aj
            aj = transport_limited()

            def combined(ac, aj):
                beta = 0.99 # smoothing factor
                # smooting the transition between Ac and Aj
                return ((ac+aj) - ((ac+aj)**2 - 4*beta*ac*aj)**0.5) / (2*beta)
            a_net = combined(ac, aj)

            #FIXME put them accordingly
            def bundle_sheath():
                alpha = 0.0001 # fraction of PSII activity in the bundle sheath cell, very low for NADP-ME types
                os = alpha * a_net / (0.047*gbs) + om # Bundle sheath O2 partial pressure, mbar
                #cbs = cm + (vp - a_net -rm) / gbs # Bundle sheath CO2 partial pressure, ubar

            return a_net

        def cost(x):
            a_net0 = x[0]
            a_net1 = c4(a_net0, press, co2, rh, tleaf)
            return (a_net0 - a_net1)**2

        #FIXME avoid passing self.stomata object to optimizer
        # iteration to obtain Cm from Ci and A, could be re-written using more efficient method like newton-raphson method
        res = scipy.optimize.minimize(cost, [0], options={'disp': True})
        a_net = res.x[0]
        self.stomata.update_stomata(pressure, co2, a_net, rh, tleaf)
        return a_net


class GasExchange:
    def __init__(self, s_type, leaf_n_content):
        self.s_type = s_type
        self.leaf_n_content = leaf_n_content

    def set_val_psil(self, pfd, tair, co2, rh, wind, press, leaf_width, leafp, et_supply):
        self.atmos = Atmosphere(pfd, tair, co2, rh, wind, press)
        self.stomata = Stomata(leaf_width)
        self.photosynthesis = Photosynthesis(self.stomata, self.leaf_n_content)

        # override GasEx() function so as to pass leaf water potential
        self._gasex_psil(leafp, et_supply)

    def _gasex_psil(self, leafp, et_supply):
        ca = self.atmos.co2
        ci = self.ci = 0.4 * ca
        self.stomata.update_boundary_layer(self.atmos.wind)

        #FIXME need initalization?
        a_net = 0.
        tleaf = self.atmos.tair
        self.stomata.update_stomata(leafp, self.atmos.co2, a_net, self.atmos.rh, tleaf)

        p = self.atmos.press
        #FIXME stomatal conductance ratio used to be 1.57, not 1.6
        a_net = (ca - ci) / self.stomata.total_resistance_co2() * p / 100.

        def pseb(pfd, press, co2, rh, leafp, tleaf):
            #FIXME minimize side-effects in _photosynthesis()
            #FIXME stomata object is the one needs to be tracked in the loop, not a_net
            a_net = self.photosynthesis.photosynthesize(pfd, press, co2, rh, leafp, tleaf)
            tleaf = self._energybalance(et_supply)
            return tleaf

        def cost(x):
            tleaf0 = x[0]
            tleaf1 = pseb(self.atmos.pfd, self.atmos.press, self.atmos.co2, self.atmos.rh, leafp, tleaf0)
            return (tleaf0 - tleaf1)**2

        res = scipy.optimize.minimize(cost, [tleaf], options={'disp': True})
        tleaf = res.x[0]

        a_net = self.photosynthesis.photosynthesize(self.atmos.pfd, self.atmos.press, self.atmos.co2, self.atmos.rh, leafp, tleaf)
        self.tleaf = tleaf

        cm = self.photosynthesis._co2_mesophyll(a_net, self.atmos.press, self.atmos.co2, self.stomata)
        self.ci = cm

        rd = self.photosynthesis._dark_respiration(tleaf)
        self.a_gross = max(0, a_net + rd) # gets negative when PFD = 0, Rd needs to be examined, 10/25/04, SK
        self.a_net = a_net

        self._evapotranspiration(self.atmos.tair, tleaf)

    def _energybalance(self, jw):
        # see Campbell and Norman (1998) pp 224-225
        # because Stefan-Boltzman constant is for unit surface area by denifition,
        # all terms including sbc are multilplied by 2 (i.e., gr, thermal radiation)
        lamda = 44000 # KJ mole-1 at 25oC
        psc = 6.66e-4
        cp = 29.3 # thermodynamic psychrometer constant and specific hear of air (J mol-1 C-1)

        #double gha, gv, gr, ghr, psc1, Ea, thermal_air, Ti, Ta;

        # variables
        ta = self.atmos.tair
        rh = self.atmos.rh
        r_abs = self.atmos.r_abs
        press = self.atmos.press

        gha = self.stomata.gb * (0.135 / 0.147) # heat conductance, gha = 1.4*.135*sqrt(u/d), u is the wind speed in m/s} Mol m-2 s-1 ?
        gv = self.stomata.total_conductance_h20()
        gr = 4 * EPS * SBC * (273 + ta)**3 / cp *2 # radiative conductance, 2 account for both sides
        ghr = gha + gr
        thermal_air = EPS * SBC * (ta + 273)**4 * 2 # emitted thermal radiation
        psc1 = psc * ghr / gv # apparent psychrometer constant

        # debug dt I commented out the changes that yang made for leaf temperature for a test. I don't think they work
        if jw == 0:
            vpd = VaporPressure.deficit(ta, rh)
            tleaf = ta + (psc1 / (VaporPressure.curve_slope(ta, press) + psc1)) * ((r_abs - thermal_air) / (ghr * cp) - vpd / (psc1 * press)) # eqn 14.6b linearized form using first order approximation of Taylor series
        else:
            tleaf = ta + (r_abs - thermal_air - lamda *jw) / (cp * ghr)
        return tleaf

    def _evapotranspiration(self, ta, tleaf):
        #variables
        ta = self.atmos.tair
        rh = self.atmos.rh
        press = self.atmos.press

        self.vpd = VaporPressure.deficit(ta, rh)

        gv = self.stomata.total_conductance_h20()
        ea = VaporPressure.ambient(ta, rh)
        es_leaf = VaporPressure.saturation(tleaf)
        et = gv * ((es_leaf - ea) / press) / (1 - (es_leaf + ea) / press)
        et = max(0., et) # 04/27/2011 dt took out the 1000 everything is moles now
        self.et = et
