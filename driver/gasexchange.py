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

    # slope of the sat vapor pressure curve: first order derivative of Es with respect to T
    @classmethod
    def curve_slope(cls, T, P):
        es = cls.saturation(T)
        b, c = cls.b, cls.c
        slope = es * (b*c)/(c+T)**2 / P
        return slope


class Atmosphere:
    def __init__(self, PFD, T_air, CO2, RH, wind, P_air):
        self.setup(PFD, T_air, CO2, RH, wind, P_air)

    def setup(self, PFD, T_air, CO2, RH, wind, P_air):
        self.PFD = PFD
        self.CO2 = CO2
        self.RH = np.clip(RH, 10, 100) / 100.
        self.T_air = T_air # C
        self.wind = wind # meters s-1
        self.P_air = P_air # kPa


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
        #FIXME initial value never used
        #self.leafp_effect = 1 # At first assume there is not drought stress, so assign 1 to leafpEffect. Yang 8/20/06

    def update_boundary_layer(self, wind):
        # maize is an amphistomatous species, assume 1:1 (adaxial:abaxial) ratio.
        sr = 1.0
        ratio = (sr + 1)**2 / (sr**2 + 1)

        # characteristic dimension of a leaf, leaf width in m
        d = self.leaf_width * 0.72

        #return 1.42 # total BLC (both sides) for LI6400 leaf chamber
        self.gb = 1.4 * 0.147 * (max(0.1, wind) / d)**0.5 * ratio
        #self.gb = (1.4 * 1.1 * 6.62 * (wind / d)**0.5 * (P_air / (R * (273.15 + T_air)))) # this is an alternative form including a multiplier for conversion from mm s-1 to mol m-2 s-1
        # 1.1 is the factor to convert from heat conductance to water vapor conductance, an avarage between still air and laminar flow (see Table 3.2, HG Jones 2014)
        # 6.62 is for laminar forced convection of air over flat plates on projected area basis
        # when all conversion is done for each surface it becomes close to 0.147 as given in Norman and Campbell
        # multiply by 1.4 for outdoor condition, Campbell and Norman (1998), p109, also see Jones 2014, pg 59 which suggest using 1.5 as this factor.
        # multiply by ratio to get the effective blc (per projected area basis), licor 6400 manual p 1-9
        return self.gb

    # stomatal conductance for water vapor in mol m-2 s-1
    def update_stomata(self, LWP, CO2, A_net, RH, T_leaf):
        # params
        g0 = self.g0
        g1 = self.g1
        gb = self.gb

        gamma = 10.0
        Cs = CO2 - (1.37 * A_net / gb) # surface CO2 in mole fraction
        if Cs <= gamma:
            Cs = gamma + 1

        m = self._leafp_effect(LWP)

        a = m * g1 * A_net / Cs
        b = g0 + gb - (m * g1 * A_net / Cs)
        c = (-RH * gb) - g0
        hs = max(np.roots([a, b, c]))
        hs = np.clip(hs, 0.3, 1.) # preventing bifurcation

        #FIXME unused?
        #es = VaporPressure.saturation(tleaf)
        #Ds = (1 - hs) * es # VPD at leaf surface
        Ds = VaporPressure.deficit(T_leaf, hs)

        gs = g0 + (g1 * m * (A_net * hs / Cs))
        gs = max(gs, g0)

        # this below is an example of how you can write temporary data to a debug window. It can be copied and
        # pasted into excel for plotting. Dennis See above where the CString object is created.
        print("gs = %f LWP = %f Ds= %f T_leaf = %f Cs = %f A_net = %f hs = %f RH = %f" % (gs, LWP, Ds, T_leaf, Cs, A_net, hs, RH))
        self.gs = gs
        return self.gs

    def _leafp_effect(self, LWP):
        # pressure - leaf water potential MPa...
        sf = 2.3 # sensitivity parameter Tuzet et al. 2003 Yang
        phyf = -1.2 # reference potential Tuzet et al. 2003 Yang
        m = (1 + np.exp(sf * phyf)) / (1 + np.exp(sf * (phyf - LWP)))
        return m

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
        self.Eac = 59400.
        self.Eao = 36000.

        self.EaVp = 75100.
        self.EaVc = 55900. # Sage (2002) JXB
        self.Eaj = 32800.

        self.Hj = 220000.
        self.Sj = 702.6

        self.Kc25 = 650. # Michaelis constant of rubisco for CO2 of C4 plants (2.5 times that of tobacco), ubar, Von Caemmerer 2000
        self.Ko25 = 450. # Michaelis constant of rubisco for O2 (2.5 times C3), mbar
        self.Kp25 = 80. # Michaelis constant for PEP caboxylase for CO2

        # Kim et al. (2007), Kim et al. (2006)
        # In von Cammerer (2000), Vpm25=120, Vcm25=60,Jm25=400
        # In Soo et al.(2006), under elevated C5O2, Vpm25=91.9, Vcm25=71.6, Jm25=354.2 YY
        self.Vpm25 = 70.
        self.Vcm25 = 50.
        self.Jm25 = 300.

        # Values in Kim (2006) are for 31C, and the values here are normalized for 25C. SK
        self.Rd25 = 2.
        self.Ear = 39800.

        #FIXME are they even used?
        #self.beta_ABA = 1.48e2 # Tardieu-Davies beta, Dewar (2002) Need the references !?
        #self.delta = -1.0
        #self.alpha_ABA = 1.0e-4
        #self.lambda_r = 4.0e-12 # Dewar's email
        #self.lambda_l = 1.0e-12
        #self.K_max = 6.67e-3 # max. xylem conductance (mol m-2 s-1 MPa-1) from root to leaf, Dewar (2002)

    def _light(self, PFD):
        #FIXME make scatt global parameter?
        scatt = 0.15 # leaf reflectance + transmittance
        f = 0.15 #spectral correction

        Ia = PFD * (1 - scatt) # absorbed irradiance
        I2 = Ia * (1 - f) / 2. # useful light absorbed by PSII
        return I2

    # mesophyll CO2 partial pressure, ubar, one may use the same value as Ci assuming infinite mesohpyle conductance
    def _co2_mesophyll(self, A_net, P_air, CO2):
        P = P_air / 100.
        Ca = CO2 * P # conversion to partial pressure
        Cm = Ca - A_net * self.stomata.total_resistance_co2() * P
        return np.clip(Cm, 0., 2*Ca)

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

    def _dark_respiration(self, T_leaf):
        return self.Rd25 * self._temperature_dependence_rate(self.Ear, T_leaf)

    def _maximum_electron_transport_rate(self, T, N):
        R = 8.314

        Tb = 25.
        K = 273.
        Tk = T + K
        Tbk = Tb + K

        Sj = self.Sj
        Hj = self.Hj

        return self.Jm25 * self._nitrogen_limited_rate(N) \
                         * self._temperature_dependence_rate(self.Eaj, T) \
                         * (1 + np.exp((Sj*Tbk - Hj) / (R*Tbk))) \
                         / (1 + np.exp((Sj*Tk  - Hj) / (R*Tk)))

    # Incident PFD, Air temp in C, CO2 in ppm, RH in percent
    #FIXME remove dependency on leaf parameters (LWP, tleaf)
    def drive(self, PFD, P_air, CO2, RH, LWP, T_leaf):
        def c4(A_net, P_air, CO2, RH, T_leaf):
            I2 = self._light(PFD)

            self.stomata.update_stomata(LWP, CO2, A_net, RH, T_leaf)
            Cm = self._co2_mesophyll(A_net, P_air, CO2)

            A_net = self.photosynthesize_c4(I2, Cm, T_leaf)
            return A_net

        def cost(x):
            A_net0 = x[0]
            A_net1 = c4(A_net0, P_air, CO2, RH, T_leaf)
            return (A_net0 - A_net1)**2

        #FIXME avoid passing self.stomata object to optimizer
        # iteration to obtain Cm from Ci and A, could be re-written using more efficient method like newton-raphson method
        res = scipy.optimize.minimize(cost, [0], options={'disp': True})
        A_net = res.x[0]
        self.stomata.update_stomata(LWP, CO2, A_net, RH, T_leaf)
        return A_net

    def photosynthesize_c4(self, I2, Cm, T_leaf):
        O = 210. # gas units are mbar
        Om = O # mesophyll O2 partial pressure

        Kp = self.Kp25 # T dependence yet to be determined
        Kc = self.Kc25 * self._temperature_dependence_rate(self.Eac, T_leaf)
        Ko = self.Ko25 * self._temperature_dependence_rate(self.Eao, T_leaf)
        Km = Kc * (1 + Om / Ko) # effective M-M constant for Kc in the presence of O2

        Rd = self._dark_respiration(T_leaf)
        Rm = 0.5 * Rd

        Vpmax = self.Vpm25 * self._nitrogen_limited_rate(self.leaf_n_content) * self._temperature_dependence_rate(self.EaVp, T_leaf)
        Vcmax = self.Vcm25 * self._nitrogen_limited_rate(self.leaf_n_content) * self._temperature_dependence_rate(self.EaVc, T_leaf)
        Jmax  = self._maximum_electron_transport_rate(T_leaf, self.leaf_n_content)

        gbs = 0.003 # bundle sheath conductance to CO2, mol m-2 s-1
        #gi = 1.0 # conductance to CO2 from intercelluar to mesophyle, mol m-2 s-1, assumed

        def enzyme_limited():
            # PEP carboxylation rate, that is the rate of C4 acid generation
            Vp1 = (Cm * Vpmax) / (Cm + Kp)
            Vp2 = Vpr = 80. # PEP regeneration limited Vp, value adopted from vC book
            Vp = max(min(Vp1, Vp2), 0.)

            #FIXME where should gamma be at?
            # half the reciprocal of rubisco specificity, to account for O2 dependence of CO2 comp point,
            # note that this become the same as that in C3 model when multiplied by [O2]
            #gamma1 = 0.193
            #gamma_star = gamma1 * Os
            #gamma = (Rd*Km + Vcmax*gamma_star) / (Vcmax - Rd)

            # Enzyme limited A (Rubisco or PEP carboxylation
            Ac1 = Vp + gbs*Cm - Rm
            #ac1 = max(0., ac1) # prevent Ac1 from being negative Yang 9/26/06
            Ac2 = Vcmax - Rd
            Ac = min(Ac1, Ac2)
            return Ac
        Ac = enzyme_limited()

        # Light and electron transport limited A mediated by J
        def transport_limited():
            theta = 0.5
            J = min(np.roots([theta, -(I2+Jmax), I2*Jmax])) # rate of electron transport
            x = 0.4 # Partitioning factor of J, yield maximal J at this value
            Aj1 = x * J/2. - Rm + gbs*Cm
            Aj2 = (1-x) * J/3. - Rd
            Aj = min(Aj1, Aj2)
            return Aj
        Aj = transport_limited()

        def combined(Ac, Aj):
            beta = 0.99 # smoothing factor
            # smooting the transition between Ac and Aj
            return ((Ac+Aj) - ((Ac+Aj)**2 - 4*beta*Ac*Aj)**0.5) / (2*beta)
        A_net = combined(Ac, Aj)

        #FIXME put them accordingly
        def bundle_sheath(A_net):
            alpha = 0.0001 # fraction of PSII activity in the bundle sheath cell, very low for NADP-ME types
            Os = alpha * A_net / (0.047*gbs) + Om # Bundle sheath O2 partial pressure, mbar
            #Cbs = Cm + (Vp - A_net - Rm) / gbs # Bundle sheath CO2 partial pressure, ubar

        return A_net


class GasExchange:
    def __init__(self, s_type, leaf_n_content):
        self.s_type = s_type
        self.leaf_n_content = leaf_n_content

    def set_val_psil(self, PFD, T_air, CO2, RH, wind, P_air, leaf_width, leafp, ET_supply):
        self.atmos = Atmosphere(PFD, T_air, CO2, RH, wind, P_air)

        stomata = Stomata(leaf_width)
        stomata.update_boundary_layer(wind)
        stomata.update_stomata(leafp, CO2, 0., RH, T_air)

        self.photosynthesis = Photosynthesis(stomata, self.leaf_n_content)

        # override GasEx() function so as to pass leaf water potential
        self._gasex_psil(self.atmos, leafp, ET_supply)

    def _gasex_psil(self, atmos, leafp, ET_supply):
        def pseb(atmos, leafp, T_leaf):
            #FIXME minimize side-effects in _photosynthesis()
            #FIXME stomata object is the one needs to be tracked in the loop, not a_net
            self.photosynthesis.drive(atmos.PFD, atmos.P_air, atmos.CO2, atmos.RH, leafp, T_leaf)
            T_leaf = self._energybalance(self.photosynthesis.stomata, atmos.T_air, atmos.RH, atmos.PFD, atmos.P_air, ET_supply)
            return T_leaf

        def cost(x):
            T_leaf0 = x[0]
            T_leaf1 = pseb(atmos, leafp, T_leaf0)
            return (T_leaf0 - T_leaf1)**2

        res = scipy.optimize.minimize(cost, [atmos.T_air], options={'disp': True})
        self.T_leaf = res.x[0]

        self.A_net = self.photosynthesis.drive(atmos.PFD, atmos.P_air, atmos.CO2, atmos.RH, leafp, self.T_leaf)
        self.Ci = self.photosynthesis._co2_mesophyll(self.A_net, atmos.P_air, atmos.CO2)

        Rd = self.photosynthesis._dark_respiration(self.T_leaf)
        self.A_gross = max(0., self.A_net + Rd) # gets negative when PFD = 0, Rd needs to be examined, 10/25/04, SK

        self._evapotranspiration(self.photosynthesis.stomata, atmos.T_air, atmos.RH, atmos.P_air, self.T_leaf)

    def _energybalance(self, stomata, T_air, RH, PFD, P_air, Jw):
        # see Campbell and Norman (1998) pp 224-225
        # because Stefan-Boltzman constant is for unit surface area by denifition,
        # all terms including sbc are multilplied by 2 (i.e., gr, thermal radiation)
        lamda = 44000 # KJ mole-1 at 25oC
        psc = 6.66e-4
        Cp = 29.3 # thermodynamic psychrometer constant and specific heat of air (J mol-1 C-1)

        epsilon = 0.97
        sbc = 5.6697e-8

        Tk = T_air + 273.

        gha = stomata.gb * (0.135 / 0.147) # heat conductance, gha = 1.4*.135*sqrt(u/d), u is the wind speed in m/s} Mol m-2 s-1 ?
        gv = stomata.total_conductance_h20()
        gr = 4 * epsilon * sbc * Tk**3 / Cp * 2 # radiative conductance, 2 account for both sides
        ghr = gha + gr
        thermal_air = epsilon * sbc * Tk**4 * 2 # emitted thermal radiation
        psc1 = psc * ghr / gv # apparent psychrometer constant

        PAR = PFD / 4.55
        # If total solar radiation unavailable, assume NIR the same energy as PAR waveband
        NIR = PAR
        scatt = 0.15
        # shortwave radiation (PAR (=0.85) + NIR (=0.15) solar radiation absorptivity of leaves: =~ 0.5
        # times 2 for projected area basis
        R_abs = (1 - scatt)*PAR + scatt*NIR + 2*(epsilon * sbc * Tk**4)

        # debug dt I commented out the changes that yang made for leaf temperature for a test. I don't think they work
        if Jw == 0:
            VPD = VaporPressure.deficit(T_air, RH)
            # eqn 14.6b linearized form using first order approximation of Taylor series
            T_leaf = T_air + (psc1 / (VaporPressure.curve_slope(T_air, P_air) + psc1)) * ((R_abs - thermal_air) / (ghr * Cp) - VPD / (psc1 * P_air))
        else:
            T_leaf = T_air + (R_abs - thermal_air - lamda * Jw) / (Cp * ghr)
        return T_leaf

    #FIXME better split into separate modules
    def _evapotranspiration(self, stomata, T_air, RH, P_air, T_leaf):
        #FIXME vpd should be in Atmosphere
        self.VPD = VaporPressure.deficit(T_air, RH)

        #FIXME et should be in Stomata
        gv = stomata.total_conductance_h20()
        ea = VaporPressure.ambient(T_air, RH)
        es_leaf = VaporPressure.saturation(T_leaf)
        ET = gv * ((es_leaf - ea) / P_air) / (1 - (es_leaf + ea) / P_air)
        ET = max(0., ET) # 04/27/2011 dt took out the 1000 everything is moles now
        self.ET = ET
