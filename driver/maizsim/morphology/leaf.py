from .organ import Organ
from ..phenology.tracker import Accumulator, BetaFunc, Q10Func, WaterStress

import numpy as np
from functools import lru_cache

class Leaf(Organ):
    def __init__(self, nodal_unit):
        super().__init__(nodal_unit.plant)
        self.nodal_unit = nodal_unit

    def setup(self):
        dt = self.p.pheno.timestep
        self._elongation_tracker = BetaFunc(R_max=1.0).use_timestep(dt)
        self._area_tracker = Accumulator().use_timestep(dt)
        self._aging_tracker = Q10Func(T_opt=self.p.pheno.optimal_temperature).use_timestep(dt)
        self._senescence_tracker = Q10Func(T_opt=self.p.pheno.optimal_temperature).use_timestep(dt)
        self._area_water_stress_tracker = WaterStress().use_timestep(dt)
        self._stay_green_water_stress_tracker = WaterStress(scale=0.5).use_timestep(dt)
        self._senescence_water_stress_tracker = WaterStress(scale=0.5).use_timestep(dt)

        # temperature at which experiments run where parameters for leaf expansion were determined
        #self._leaf_calibrated_temperature = self.p.pheno.calibrated_tempreature

        #FIXME the original code actually doesn't do this, always use the generic leaf number
        # potential growth is fixed by the number of total leaves at initiation
        #self._total_leaves_at_initiation = self.p.pheno.leaves_total

    #############
    # Constants #
    #############

    @property
    def rank(self):
        return self.nodal_unit.rank

    # cm dd-1 Fournier and Andrieu 1998 Pg239.
    # This is the "potential" elongation rate with no water stress Yang
    @property
    def elongation_rate(self):
        return 0.564

    # max elongation rate (cm per day) at optipmal temperature
    # (Topt: 31C with Tbase = 9.8C using 0.564 cm/dd rate from Fournier 1998 paper above
    @property
    def maximum_elongation_rate(self):
        return 12.0

    # leaf lamina width to length ratio
    @property
    def width_to_length_ratio(self):
        return 0.106

    # leaf area coeff with respect to L*W (A_LW)
    @property
    def area_ratio(self):
        return 0.75

    # staygreen trait of the hybrid
    # stay green for this value times growth period after peaking before senescence begins
    # An analogy for this is that with no other stresses involved,
    # it takes 15 years to grow up, stays active for 60 years,
    # and age the last 15 year if it were for a 90 year life span creature.
    # Once fully grown, the clock works differently so that the hotter it is quicker it ages
    @property
    def stay_green(self):
        return 3.5

    #############
    # Variables #
    #############

    @property
    @lru_cache()
    def _extra_leaves(self):
        #FIXME used to leaves_total, not potential
        return self.p.pheno.leaves_potential - self.p.pheno.leaves_generic

    @lru_cache()
    def _maximum_length(self, extra_leaves):
        LM_min = 115
        k = 24.0
        return np.sqrt(LM_min**2 + k * extra_leaves)

    @property
    @lru_cache()
    def maximum_length(self):
        return self._maximum_length(self._extra_leaves)

    @property
    @lru_cache()
    def maximum_width(self):
        # Fournier and Andrieu(1998) Pg242 YY
        return self.maximum_length * self.width_to_length_ratio

    #TODO better name, shared by growth_duration and pontential_area
    @lru_cache()
    def _rank_effect(self, rank, leaves, weight):
        n_m = 5.93 + 0.33 * leaves # the rank of the largest leaf. YY
        a = (-10.61 + 0.25 * leaves) * weight
        b = (-5.99 + 0.27 * leaves) * weight
        # equation 7 in Fournier and Andrieu (1998). YY

        # equa 8(b)(Actually eqn 6? - eqn 8 deals with leaf age - DT)
        # in Fournier and Andrieu(1998). YY
        scale = rank / n_m - 1
        return np.exp(a * scale**2 + b * scale**3)

    def rank_effect(self, weight=1):
        #TODO should be a plant parameter not leaf (?)
        #FIXME leaves_potential used to be leaves_total
        return self._rank_effect(self.rank, self.p.pheno.leaves_potential, weight)

    @property
    @lru_cache()
    def potential_length(self):
        return self.maximum_length * self.rank_effect(weight=0.5)

    # from CLeaf::calc_dimensions()
    # LM_min is a length characteristic of the longest leaf,in Fournier and Andrieu 1998, it was 90 cm
    # LA_max is a fn of leaf no (Birch et al, 1998 fig 4) with largest reported value near 1000cm2. This is implemented as lfno_effect below, SK
    # LM_min of 115cm gives LA of largest leaf 1050cm2 when totalLeaves are 25 and Nt=Ng, SK 1-20-12
    # Without lfno_effect, it can be set to 97cm for the largest leaf area to be at 750 cm2 with Nt ~= Ng (Lmax*Wmax*0.75) based on Muchow, Sinclair, & Bennet (1990), SK 1-18-2012
    # Eventually, this needs to be a cultivar parameter and included in input file, SK 1-18-12
    # the unit of k is cm^2 (Fournier and Andrieu 1998 Pg239). YY
    # L_max is the length of the largest leaf when grown at T_peak. Here we assume LM_min is determined at growing Topt with minmal (generic) leaf no, SK 8/2011
    # If this routine runs before TI, totalLeaves = genericLeafNo, and needs to be run with each update until TI and total leaves are finalized, SK
    @property
    @lru_cache()
    def growth_duration(self):
        # shortest possible linear phase duration in physiological time (days instead of GDD) modified
        return self.potential_length / self.maximum_elongation_rate

    @property
    @lru_cache()
    def phase1_delay(self):
        # not used in MAIZSIM because LTAR is used to initiate leaf growth.
        # Fournier's value : -5.16+1.94*rank;equa 11 Fournier and Andrieu(1998) YY, This is in plastochron unit
        return np.fmax(0, -5.16 + 1.94 * self.rank)

    @lru_cache()
    def _leaf_number_effect(self, leaves):
        # Fig 4 of Birch et al. (1998)
        return np.clip(np.exp(-1.17 + 0.047 * leaves), 0.5, 1.0)

    @property
    @lru_cache()
    def leaf_number_effect(self):
        #FIXME used to leaves_total, not potential
        return self._leaf_number_effect(self.p.pheno.leaves_potential)

    @property
    @lru_cache()
    def potential_area(self):
        # daughtry and hollinger (1984) Fournier and Andrieu(1998) Pg242 YY
        maximum_area = self.maximum_length * self.maximum_width * self.area_ratio

        # equa 6. Fournier and Andrieu(1998) multiplied by Birch et al. (1998) leaf no effect
        # LA_max the area of the largest leaf
        # PotentialArea potential final area of a leaf with rank "n". YY
        return maximum_area * self.leaf_number_effect * self.rank_effect(weight=1)

    @property
    def green_ratio(self):
        return 1 - self.senescence_ratio

    @property
    def green_area(self):
        return self.green_ratio * self.area

    @property
    def elongation_age(self):
        #TODO implement Parent and Tardieu (2011, 2012) approach for leaf elongation in response to T and VPD, and normalized at 20C, SK, Nov 2012
        # elongAge indicates where it is now along the elongation stage or duration.
        # duration is determined by totallengh/maxElongRate which gives the shortest duration to reach full elongation in the unit of days.
        #FIXME no need to check here, as it will be compared against duration later anyways
        #return np.fmin(self._elongation_tracker.rate, self.growth_duration)
        return self._elongation_tracker.rate

    @lru_cache()
    def _temperature_effect(self, T_grow, T_peak=18.7, T_base=8.0):
        # T_peak is the optimal growth temperature at which the potential leaf size determined in calc_mophology achieved.
        # Similar concept to fig 3 of Fournier and Andreiu (1998)

        # phyllochron corresponds to PHY in Lizaso (2003)
        # phyllochron needed for next leaf appearance in degree days (GDD8) - 08/16/11, SK.
        #phyllochron = (dv->get_T_Opt()- Tb)/(dv->get_Rmax_LTAR());

        T_ratio = (T_grow - T_base) / (T_peak - T_base)
        # final leaf size is adjusted by growth temperature determining cell size during elongation
        return np.fmax(0, T_ratio * np.exp(1 - T_ratio))

    @property
    def temperature_effect(self):
        return self._temperature_effect(T_grow=self.p.pheno.growing_temperature)

    #TODO confirm if it really means the elongation rate
    @property
    def elongation_rate(self):
        t = self.elongation_age
        t_e = self.growth_duration
        t = np.fmin(t, t_e)
        #HACK can be more simplified
        #a = 1.5 / t_e
        #b = 4 * t * (t_e - t) / t_e**2
        #return a * b
        t_m = t_e / 2
        a = (2*t_e - t_m) / (t_e * (t_e - t_m)) * (t_m / t_e)**(t_m / (t_e - t_m))
        b = np.fmax(0, (t_e - t) / (t_e - t_m) * (t / t_m)**(t_m / (t_e - t_m)))
        return np.fmax(0, a * b)

    @property
    def potential_area_increase(self):
        # time step as day fraction
        timestep = self._elongation_tracker.timestep

        ##area = np.fmax(0, water_effect * T_effect * self.potential_area * (1 + (t_e - self.elongation_age) / (t_e - t_m)) * (self.elongation_age / t_e)**(t_e / (t_e - t_m)))
        #maximum_expansion_rate = T_effect * self.potential_area * (2*t_e - t_m) / (t_e * (t_e - t_m)) * (t_m / t_e)**(t_m / (t_e - t_m))
        # potential leaf area increase without any limitations
        #return np.fmax(0, maximum_expansion_rate * np.fmax(0, (t_e - self.elongation_age) / (t_e - t_m) * (self.elongation_age / t_m)**(t_m / (t_e - t_m))) * timestep)
        return self.temperature_effect * self.elongation_rate * timestep * self.potential_area

    # create a function which simulates the reducing in leaf expansion rate
    # when predawn leaf water potential decreases. Parameterization of rf_psil
    # and rf_sensitivity are done with the data from Boyer (1970) and Tanguilig et al (1987) YY
    @lru_cache()
    def _water_potential_effect(self, psi_predawn, threshold):
        #psi_predawn = self.p.soil.WP_leaf_predawn
        psi_th = threshold # threshold wp below which stress effect shows up

        # DT Oct 10, 2012 changed this so it was not as sensitive to stress near -0.5 lwp
        # SK Sept 16, 2014 recalibrated/rescaled parameter estimates in Yang's paper. The scale of Boyer data wasn't set correctly
        # sensitivity = 1.92, LeafWPhalf = -1.86, the sensitivity parameter may be raised by 0.3 to 0.5 to make it less sensitivy at high LWP, SK
        s_f = 0.4258 # 0.5
        psi_f = -1.4251 # -1.0
        return np.fmin(1.0, (1 + np.exp(s_f * psi_f)) / (1 + np.exp(s_f * (psi_f - (psi_predawn - psi_th)))))

    def water_potential_effect(self, threshold):
        return self._water_potential_effect(self.p.soil.WP_leaf_predawn, threshold)

    @property
    def actual_area_increase(self):
        # See Kim et al. (2012) Agro J. for more information on how this relationship has been derermined basned on multiple studies and is applicable across environments
        water_effect = self.water_potential_effect(-0.8657)

        # place holder
        carbon_effect = 1.0

        # growth temperature effect is included in determining potential area
        return np.fmin(water_effect, carbon_effect) * self.potential_area_increase

    @property
    def relative_area_increase(self):
        # adapted from CPlant::calcPerLeafRelativeAreaIncrease()
        return self.potential_area_increase / self.nodal_unit.plant.potential_leaf_area_increase

    @property
    # actual area
    def area(self):
        return self._area_tracker.rate

    @property
    def stay_green_water_stress_duration(self):
        return self._stay_green_water_stress_tracker.rate

    @property
    def stay_green_duration(self):
        # SK 8/20/10: as in Sinclair and Horie, 1989 Crop sciences, N availability index scaled between 0 and 1 based on
        #nitrogen_index = np.fmax(0, (2 / (1 + np.exp(-2.9 * (self.g_content - 0.25))) - 1))
        return np.fmax(0, self.stay_green * self.growth_duration - self.stay_green_water_stress_duration)

    @property
    def active_age(self):
        # Assumes physiological time for senescence is the same as that for growth though this may be adjusted by stayGreen trait
        # a peaked fn like beta fn not used here because aging should accelerate with increasing T not slowing down at very high T like growth,
        # instead a q10 fn normalized to be 1 at T_opt is used, this means above Top aging accelerates.
        #FIXME no need to check here, as it will be compared against duration later anyways
        #return np.fmin(self._aging_tracker.rate, self.stay_green_duration)
        return self._aging_tracker.rate

    @property
    def senescence_water_stress_duration(self):
        return self._senescence_water_stress_tracker.rate

    @property
    def senescence_duration(self):
        # end of growth period, time to maturity
        return np.fmax(0, self.growth_duration - self.senescence_water_stress_duration)

    @property
    def senescence_age(self):
        #FIXME no need to check here, as it will be compared against duration later anyways
        #return np.fmin(self._senescence_tracker.rate, self.senescence_duration)
        return self._senescence_tracker.rate

    @property
    #TODO confirm if it really means the senescence ratio, not rate
    def senescence_ratio(self):
        t = self.senescence_age
        t_e = self.senescence_duration
        if t >= t_e:
            return 1
        else:
            t_m = t_e / 2
            r = (1 + (t_e - t) / (t_e - t_m)) * (t / t_e)**(t_e / (t_e - t_m))
            return np.clip(r, 0, 1)

    @property
    def senescent_area(self):
        # Leaf senescence accelerates with drought and heat. see http://www.agry.purdue.edu/ext/corn/news/timeless/TopLeafDeath.html
        #timestep = self._elongation_tracker.timestep
        #rate = self._growth_rate(self.senescence_age, self.senescence_duration)
        #return rate * timestep * self.area
        return self.senescence_ratio * self.area

    @property
    def specific_leaf_area(self):
        # temporary for now - it should vary by age. Value comes from some of Soo's work
        #return 200.0
        return self.area / self.mass


    # Nitrogen

    @property
    def nitrogen(self):
        #TODO is this default value needed?
        # no N stress
        #return 3.0
        return self.p.nitrogen.leaf_content

    ##########
    # States #
    ##########

    @property
    def initiated(self):
        # no explicit initialize() here
        return True

    @property
    def appeared(self):
        return self.rank <= self.p.pheno.leaves_appeared

    @property
    def growing(self):
        return self.appeared and not self.mature

    @property
    def mature(self):
        return self.elongation_age >= self.growth_duration or self.area >= self.potential_area

    @property
    def aging(self):
        return self.active_age >= self.stay_green_duration

    @property
    def dead(self):
        #return self.senescent_area >= self.area
        return self.senescence_ratio >= 1
        #return self.senescence_age >= self.senescence_duration?

    @property
    def dropped(self):
        return self.mature and self.dead

    ##########
    # Update #
    ##########

    #FIXME signature mismatch with Organ: T vs predawn_lwp
    def update(self):
        super().update()
        self.expand()
        self.senescence()
        self.water_stress()

    # leaf expansiopn rate based on a determinate sigmoid function by Yin et al. (2003)
    def expand(self):
        if self.appeared and not self.mature:
            self._elongation_tracker.update(self.p.pheno.temperature)
            self._area_tracker.update(self.actual_area_increase)

    def senescence(self):
        T = self.p.pheno.temperature
        if not self.aging:
            self._aging_tracker.update(T)
        elif not self.dead:
            self._senescence_tracker.update(T)

    def water_stress(self):
        if self.mature:
            # One day of cumulative severe water stress (i.e., water_effect = 0.0 around -4MPa) would result in a reduction of leaf lifespan in relation staygreeness and growthDuration, SK
            # if scale is 1.0, one day of severe water stress shortens one day of stayGreenDuration
            self._stay_green_water_stress_tracker.update(self.water_potential_effect(-4.0))
        if self.aging:
            # if scale is 0.5, one day of severe water stress at predawn shortens one half day of agingDuration
            self._senescence_water_stress_tracker.update(self.water_potential_effect(-4.0))
        #TODO handle dead leaf?
        #elif self.dead:
        #TODO handle after aging?
