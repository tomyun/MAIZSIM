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

R = 8.314
MAXITER = 200
#FIXME are they parameters?
EPS = 0.97
SBC = 5.6697e-8

class Atmosphere:
    def __init__(self):
        pass

    @classmethod
    def saturation_vapor_pressure(cls, t):
        # Campbell and Norman (1998), p 41 Saturation vapor pressure in kPa
        a = 0.611 # kPa
        b = 17.502
        c = 240.97 # C

        #FIXME August-Roche-Magnus formula gives slightly different parameters
        # https://en.wikipedia.org/wiki/Clausius–Clapeyron_relation
        #a = 0.61094 # kPa
        #b = 17.625
        #c = 243.04 # C

        return a*np.exp((b*t)/(c+t))

    @classmethod
    def vapor_pressure(cls, t, rh):
        es = cls.saturation_vapor_pressure(t)
        return es * rh

    @classmethod
    def vapor_pressure_deficit(cls, t, rh):
        es = cls.saturation_vapor_pressure(t)
        return es * (1 - rh)

class Stomata:
    def __init__(self):
        self.setup()

    def setup(self):
        # in P. J. Sellers, et al.Science 275, 502 (1997)
        # g0 is b, of which the value for c4 plant is 0.04
        # and g1 is m, of which the value for c4 plant is about 4 YY
        self.g0 = 0.04
        self.g1 = 4.0

    @classmethod
    def boundary_layer_conductance(cls, width, wind):
        # maize is an amphistomatous species, assume 1:1 (adaxial:abaxial) ratio.
        sr = 1.0
        ratio = (sr + 1)**2 / (sr**2 + 1)

        # characteristic dimension of a leaf, leaf width in m
        d = width * 0.72

        #return 1.42 # total BLC (both sides) for LI6400 leaf chamber
        return 1.4 * 0.147 * (max(0.1, wind) / d)**0.5 * ratio
        # return (1.4 * 1.1 * 6.62 * (wind / d)**0.5 * (press / (R * (273.15 + tair)))) # this is an alternative form including a multiplier for conversion from mm s-1 to mol m-2 s-1
        # 1.1 is the factor to convert from heat conductance to water vapor conductance, an avarage between still air and laminar flow (see Table 3.2, HG Jones 2014)
        # 6.62 is for laminar forced convection of air over flat plates on projected area basis
        # when all conversion is done for each surface it becomes close to 0.147 as given in Norman and Campbell
        # multiply by 1.4 for outdoor condition, Campbell and Norman (1998), p109, also see Jones 2014, pg 59 which suggest using 1.5 as this factor.
        # multiply by ratio to get the effective blc (per projected area basis), licor 6400 manual p 1-9

    # stomatal conductance for water vapor in mol m-2 s-1
    def stomatal_conductance(self, pressure, co2, a_net, gb, rh, tleaf):
        # params
        g0 = self.g0
        g1 = self.g1

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

        #es = Atmosphere.saturation_vapor_pressure(tleaf)
        #ds = (1 - hs) * es # VPD at leaf surface
        ds = Atmosphere.vapor_pressure_deficit(tleaf, hs)
        tmp = g0 + (g1 *temp * (a_net * hs / cs))
        tmp = max(tmp, g0)

        # this below is an example of how you can write temporary data to a debug window. It can be copied and
        # pasted into excel for plotting. Dennis See above where the CString object is created.
        print("tmp = %f pressure = %f Ds= %f Tleaf = %f Cs = %f Anet = %f hs = %f RH = %f" % (tmp, pressure, ds, tleaf, cs, a_net, hs, rh))
        return tmp

    def _set_leafp_effect(self, pressure):
        # pressure - leaf water potential MPa...
        sf = 2.3 # sensitivity parameter Tuzet et al. 2003 Yang
        phyf = -1.2 # reference potential Tuzet et al. 2003 Yang
        self.leafp_effect = (1 + np.exp(sf * phyf)) / (1 + np.exp(sf * (phyf - pressure)))
        return self.leafp_effect

    @classmethod
    def total_conductance(cls, gb, gs):
        return gs * gb / (gs + gb)


class GasExchange:
    def __init__(self, s_type, n_content):
        self.s_type = s_type
        self.leaf_n_content = n_content

        self.stomata = Stomata()

    def set_val_psil(self, pfd, tair, co2, rh, wind, press, width, leafp, et_supply):
        scatt = 0.15
        self.pfd = pfd
        par = pfd / 4.55
        # If total solar radiation unavailable, assume NIR the same energy as PAR waveband
        nir = par
        # times 2 for projected area basis
        self.r_abs = (1 - scatt)*par + 0.15*nir + 2*(EPS * SBC * (tair+273)**4)
        # shortwave radiation (PAR (=0.85) + NIR (=0.15) solar radiation absorptivity of leaves: =~ 0.5
        self.co2 = co2
        self.rh = np.clip(rh, 10, 100) / 100.
        self.tair = tair # C
        self.width = width / 100. # meters
        self.wind = wind # meters s-1
        self.press = press # kPa
        self.leafp = leafp
        self.leafp_effect = 1 # At first assume there is not drought stress, so assign 1 to leafpEffect. Yang 8/20/06

        self._setup_parms()

        # override GasEx() function so as to pass leaf water potential
        self._gasex_psil(leafp, et_supply)

    def _setup_parms(self):
        self.eavp = 75100
        self.eavc = 55900 # Sage (2002) JXB
        self.eaj = 32800
        self.hj = 220000
        self.sj = 702.6

        # Kim et al. (2007), Kim et al. (2006)
        # In von Cammerer (2000), Vpm25=120, Vcm25=60,Jm25=400
        # In Soo et al.(2006), under elevated C5O2, Vpm25=91.9, Vcm25=71.6, Jm25=354.2 YY
        self.vpm25 = 70.0
        self.vcm25 = 50.0
        self.jm25 = 300.0

        # Values in Kim (2006) are for 31C, and the values here are normalized for 25C. SK
        self.rd25 = 2.0
        self.ear = 39800

        self.beta_aba = 1.48e2 # Tardieu-Davies beta, Dewar (2002) Need the references !?
        self.delta = -1.0
        self.alpha_aba = 1.0e-4
        self.lambda_r = 4.0e-12 # Dewar's email
        self.lambda_l = 1.0e-12
        self.k_max = 6.67e-3 # max. xylem conductance (mol m-2 s-1 MPa-1) from root to leaf, Dewar (2002)

    def _gasex_psil(self, leafp, et_supply):
        ca = self.co2
        ci = self.ci = 0.4 * ca
        gb = self.gb = self.stomata.boundary_layer_conductance(self.width, self.wind)

        #FIXME need initalization?
        self.a_net = 0.
        self.tleaf = self.tair
        gs = self.gs = self.stomata.stomatal_conductance(leafp, self.co2, self.a_net, gb, self.rh, self.tleaf)

        p = self.press
        self.a_net = (ca - ci) / (1.57 / gs + 1.37 / gb) * p / 100.

        i = 1
        tleaf_old = 0.
        while abs(tleaf_old - self.tleaf) > 0.01 and i < MAXITER:
            tleaf_old = self.tleaf
            self._photosynthesis(leafp, self.tleaf)
            self.tleaf = self._energybalance(et_supply)
            self._evapotranspiration(self.tair, self.tleaf)
            i += 1
            #FIXME remove
            self.iter2 = i

    # Incident PFD, Air temp in C, CO2 in ppm, RH in percent
    def _photosynthesis(self, pressure, tleaf):
        f = 0.15 #spectral correction
        o = 210 # gas units are mbar
        theta = 0.5
        #FIXME make scatt global parameter?
        scatt = 0.15 # leaf reflectance + transmittance
        kc25 = 650 # Michaelis constant of rubisco for CO2 of C4 plants (2.5 times that of tobacco), ubar, Von Caemmerer 2000
        ko25 = 450 # Michaelis constant of rubisco for O2 (2.5 times C3), mbar
        kp25 = 80 # Michaelis constant for PEP caboxylase for CO2
        eac = 59400
        eao = 36000 # activation energy values
        vpr = 80 # PEP regeneration limited Vp, value adopted from vC book
        gbs = 0.003 # bundle sheath conductance to CO2, mol m-2 s-1
        x = 0.4 # Partitioning factor of J, yield maximal J at this value
        alpha = 0.0001 # fraction of PSII activity in the bundle sheath cell, very low for NADP-ME types
        gi = 1.0 # conductance to CO2 from intercelluar to mesophyle, mol m-2 s-1, assumed
        beta = 0.99 # smoothing factor
        gamma1 = 0.193 # half the reciprocal of rubisco specificity, to account for O2 dependence of CO2 comp point, note that this become the same as that in C3 model when multiplied by [O2]

        #FIXME is it necessary?
        # Reset values changed for N status
        self._setup_parms()

        # variables
        pfd = self.pfd
        #tleaf = self.tleaf
        press = self.press
        co2 = self.co2

        # params
        rd25 = self.rd25
        ear = self.ear
        vcm25 = self.vcm25
        jm25 = self.jm25
        vpm25 = self.vpm25
        eavp = self.eavp
        eavc = self.eavc
        eaj = self.eaj
        sj = self.sj
        hj = self.hj

        #double Kp, Kc, Ko, Km, Ia, I2, Tk, Vpmax, Jmax, Vcmax, Om, Rm, J, Ac1, Ac2, Ac, Aj1,
        #       Aj2, Aj, Vp1, Vp2, Vp, P,  Ca, Cm, Cm_last, Cm_next,
        #       Os, GammaStar, Gamma;

        # Light response function parameters
        ia = pfd * (1 - scatt) # absorbed irradiance
        i2 = ia * (1 - f) / 2. # useful light absorbed by PSII

        # other input parameters and constants
        tk = tleaf + 273.0
        p = press / 100.

        #FIXME name it properly
        tleaf_ratio = (tleaf - 25) / (298 * R * (tleaf + 273))

        ca = co2 * p # conversion to partial pressure
        om = o # mesophyle O2 partial pressure
        kp = kp25 # T dependence yet to be determined
        kc = kc25 * np.exp(eac * tleaf_ratio)
        ko = ko25 * np.exp(eao * tleaf_ratio)
        km = kc * (1 + om / ko) # effective M-M constant for Kc in the presence of O2
        rd = rd25 * np.exp(ear * tleaf_ratio)

        critical_nitrogen = max(0.25, self.leaf_n_content)
        cnr = 2 / (1 + np.exp(-2.9 * (critical_nitrogen - 0.25))) - 1
        vcm25 *= cnr
        jm25  *= cnr
        vpm25 *= cnr
        # in Sinclair and Horie, 1989 Crop sciences, it is 4 and 0.2
        # In J Vos. et al. Field Crop study, 2005, it is 2.9 and 0.25
        # In Lindquist, weed science, 2001, it is 3.689 and 0.5

        vpmax = vpm25 * np.exp(eavp * tleaf_ratio)
        vcmax = vcm25 * np.exp(eavc * tleaf_ratio)
        jmax  = jm25  * np.exp(((tk - 298) * eaj) / (R * tk * 298)) \
                * (1 + np.exp((sj * 298 - hj) / (R * 298))) \
                / (1 + np.exp((sj * tk - hj) / (R * tk)))
        rm = 0.5 * rd

        # mesophyll CO2 partial pressure, ubar, one may use the same value as Ci assuming infinite mesohpyle conductance
        def co2_mesophyll(ca, a_net, gs, gb):
            cm = ca - a_net * (1.6 / gs + 1.37 / gb) * p
            return np.clip(cm, 0., 2*ca)

        def c4(a_net, gb, co2, rh, tleaf):
            gs = self.stomata.stomatal_conductance(pressure, co2, a_net, gb, rh, tleaf)
            cm = co2_mesophyll(ca, a_net, gs, gb)

            # PEP carboxylation rate, that is the rate of C4 acid generation
            vp1 = (cm * vpmax) / (cm + kp)
            vp2 = vpr
            vp = max(min(vp1, vp2), 0)

            # Enzyme limited A (Rubisco or PEP carboxylation
            ac1 = vp + gbs*cm - rm
            #ac1 = max(0, ac1) # prevent Ac1 from being negative Yang 9/26/06
            ac2 = vcmax - rd
            ac = min(ac1, ac2)

            # Light and electron transport limited A mediated by J
            #j = QuadSolnLower(theta, -(i2+jmax), i2*jmax) # rate of electron transport
            j = min(np.roots([theta, -(i2+jmax), i2*jmax])) # rate of electron transport
            aj1 = x*j/2. - rm + gbs*cm
            aj2 = (1 - x)*j/3. - rd
            aj = min(aj1, aj2)

            return ((ac+aj) - ((ac+aj)**2 - 4*beta*ac*aj)**0.5) / (2*beta) # smooting the transition between Ac and Aj

        def cost(x):
            a_net0 = x[0]
            a_net1 = c4(a_net0, self.gb, self.co2, self.rh, self.tleaf)
            return (a_net0 - a_net1)**2

        # iteration to obtain Cm from Ci and A, could be re-written using more efficient method like newton-raphson method
        res = scipy.optimize.minimize(cost, [self.a_net], options={'disp': True})
        self.a_net = res.x[0]
        self.gs = self.stomata.stomatal_conductance(pressure, self.co2, self.a_net, self.gb, self.rh, self.tleaf)
        cm = co2_mesophyll(ca, self.a_net, self.gs, self.gb)

        os = alpha * self.a_net / (0.047*gbs) + om # Bundle sheath O2 partial pressure, mbar
        #cbs = cm + (vp - a_net -rm) / gbs # Bundle sheath CO2 partial pressure, ubar
        gamma_star = gamma1 * os
        gamma = (rd*km + vcmax*gamma_star) / (vcmax - rd)
        self.ci = cm
        self.a_gross = max(0, self.a_net + rd) # gets negative when PFD = 0, Rd needs to be examined, 10/25/04, SK

    def _energybalance(self, jw):
        # see Campbell and Norman (1998) pp 224-225
        # because Stefan-Boltzman constant is for unit surface area by denifition,
        # all terms including sbc are multilplied by 2 (i.e., gr, thermal radiation)
        lamda = 44000 # KJ mole-1 at 25oC
        psc = 6.66e-4
        cp = 29.3 # thermodynamic psychrometer constant and specific hear of air (J mol-1 C-1)

        #double gha, gv, gr, ghr, psc1, Ea, thermal_air, Ti, Ta;

        # variables
        ta = self.tair
        #ti = self.tleaf
        gb = self.gb
        gs = self.gs
        rh = self.rh
        r_abs = self.r_abs
        press = self.press

        gha = gb * (0.135 / 0.147) # heat conductance, gha = 1.4*.135*sqrt(u/d), u is the wind speed in m/s} Mol m-2 s-1 ?
        gv = self.stomata.total_conductance(gb, gs)
        gr = 4 * EPS * SBC * (273 + ta)**3 / cp *2 # radiative conductance, 2 account for both sides
        ghr = gha + gr
        thermal_air = EPS * SBC * (ta + 273)**4 * 2 # emitted thermal radiation
        psc1 = psc * ghr / gv # apparent psychrometer constant

        # debug dt I commented out the changes that yang made for leaf temperature for a test. I don't think they work
        if jw == 0:
            vpd = Atmosphere.vapor_pressure_deficit(ta, rh)
            tleaf = ta + (psc1 / (self._slope(ta) + psc1)) * ((r_abs - thermal_air) / (ghr * cp) - vpd / (psc1 * press)) # eqn 14.6b linearized form using first order approximation of Taylor series
        else:
            tleaf = ta + (r_abs - thermal_air - lamda *jw) / (cp * ghr)
        return tleaf

    def _evapotranspiration(self, ta, tleaf):
        #variables
        ta = self.tair
        gb = self.gb
        gs = self.gs
        rh = self.rh
        press = self.press

        self.vpd = Atmosphere.vapor_pressure_deficit(ta, rh)

        gv = self.stomata.total_conductance(gb, gs)
        ea = Atmosphere.vapor_pressure(ta, rh) # ambient vapor pressure
        es_leaf = Atmosphere.saturation_vapor_pressure(tleaf)
        et = gv * ((es_leaf - ea) / press) / (1 - (es_leaf + ea) / press)
        et = max(0., et) # 04/27/2011 dt took out the 1000 everything is moles now
        self.et = et

    # slope of the sat vapor pressure curve: first order derivative of Es with respect to T
    def _slope(self, t):
        # variables
        press = self.press

        # units of b and c are  degrees C
        b = 17.502
        c = 240.97

        es = Atmosphere.saturation_vapor_pressure(t)
        slope = es * (b*c) / (c + t)**2 / press
        return slope
