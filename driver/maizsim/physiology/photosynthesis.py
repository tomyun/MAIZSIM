from .trait import Trait
from .gasexchange import GasExchange
from ..atmosphere import Sun
from ..morphology import Radiation
from ..morphology.weight import Weight

import numpy as np

#TODO rename to CarbonAssimilation or so? could be consistently named as CarbonPartition, CarbonAllocation...
class Photosynthesis(Trait):
    def setup(self):
        self.sunlit = GasExchange('Sunlit')
        self.shaded = GasExchange('Shaded')
        self.radiation = None

    # calc_gas_exchange() from Plant
    def update(self):
        #tau = 0.50 # atmospheric transmittance, to be implemented as a variable => done

        LAF = 1.37 # leaf angle factor for corn leaves, Campbell and Norman (1998)
        leaf_width = 5.0 # to be calculated when implemented for individal leaves
        LAI = self.p.area.leaf_area_index

        #TODO how do we get LeafWP and ET_supply?
        LWP = self.p.soil.WP_leaf
        ET_supply = self.p.water.supply * self.p.initials.plant_density / 3600 / 18.01 / LAI

        # jday = Timer.julian_day_from_datetime(self.p.weather.time)
        # jhour = Timer.julian_hour_from_datetime(self.p.weather.time)

        #TODO integrate lightenv with Atmosphere class?
        #TODO lightenv.dll needs to be translated to C++. It slows down the execution, 3/16/05, SK
        # self.lightenv.radTrans2(
        #     jday, jhour,
        #     self.p.initials.latitude, self.p.initials.longitude,
        #     self.p.weather.sol_rad, self.p.weather.PFD,
        #     LAI, LAF
        # )
        # temp7 = lightenv.getNIRtot()

        sun = Sun(self.p.weather.time, self.p.initials.latitude, self.p.initials.longitude, PAR=self.p.weather.PFD)
        self.radiation = Radiation(sun, leaf_area_index=LAI, leaf_angle_factor=LAF)

        # Calculating transpiration and photosynthesis without stomatal control Y
        # call SetVal_NC()

        #import pdb; pdb.set_trace()

        # Calculating transpiration and photosynthesis with stomatal controlled by leaf water potential LeafWP Y
        sunlit_weather = self.p.weather.copy()
        sunlit_weather.PFD = self.radiation.irradiance_Q_sunlit()
        self.sunlit.setup(sunlit_weather, self.p.soil, self.p.nitrogen.leaf_content, leaf_width, ET_supply)

        shaded_weather = self.p.weather.copy()
        shaded_weather.PFD = self.radiation.irradiance_Q_shaded()
        self.shaded.setup(shaded_weather, self.p.soil, self.p.nitrogen.leaf_content, leaf_width, ET_supply)

    @property
    def sunlit_leaf_area_index(self):
        return self.radiation.sunlit_leaf_area_index if self.radiation else 0

    @property
    def shaded_leaf_area_index(self):
        return self.radiation.shaded_leaf_area_index if self.radiation else 0

    @property
    def leaf_area_index_array(self):
        return np.array([
            self.sunlit_leaf_area_index,
            self.shaded_leaf_area_index,
        ])

    def _weighted(self, array):
        return self.leaf_area_index_array.dot(array)

    @property
    def sunlit_irradiance(self):
        return self.radiation.irradiance_Q_sunlit() if self.radiation else 0

    @property
    def shaded_irradiance(self):
        return self.radiation.irradiance_Q_shaded() if self.radiation else 0

    @property
    def gross_array(self):
        return np.array([
            self.sunlit.A_gross,
            self.shaded.A_gross,
        ])

    @property
    def net_array(self):
        return np.array([
            self.sunlit.A_net,
            self.shaded.A_net,
        ])

    @property
    def evapotranspiration_array(self):
        return np.array([
            self.sunlit.ET,
            self.shaded.ET,
        ])

    @property
    def temperature_array(self):
        return np.array([
            self.sunlit.T_leaf,
            self.shaded.T_leaf,
        ])

    @property
    def conductance_array(self):
        return np.array([
            self.sunlit.gs,
            self.shaded.gs,
        ])

    @property
    def gross_CO2_umol_per_m2_s(self):
        return self._weighted(self.gross_array)

    # plantsPerMeterSquare units are umol CO2 m-2 ground s-1
    # in the following we convert to g C plant-1 per hour
    # photosynthesis_gross is umol CO2 m-2 leaf s-1

    @property
    def net_CO2_umol_per_m2_s(self):
        # grams CO2 per plant per hour
        return self._weighted(self.net_array)

    @property
    def transpiration_H2O_mol_per_m2_s(self):
        #TODO need to save this?
        # when outputting the previous step transpiration is compared to the current step's water uptake
        #self.transpiration_old = self.transpiration
        #FIXME need to check if LAIs are negative?
        #transpiration = sunlit.ET * max(0, sunlit_LAI) + shaded.ET * max(0, shaded_LAI)
        return self._weighted(self.evapotranspiration_array)

    #TODO consolidate unit conversions somewhere else

    @property
    def _mol_per_umol(self):
        return 1 / 1e6

    @property
    def _plant_per_m2(self):
        return 1 / self.p.initials.plant_density

    @property
    def _min_step_per_sec(self):
        return 60 * self.p.initials.timestep

    # final values

    @property
    def assimilation(self):
        # grams CO2 per plant per hour
        return np.prod([
            self.gross_CO2_umol_per_m2_s,
            self._mol_per_umol,
            self._plant_per_m2,
            self._min_step_per_sec,
            Weight.CO2,
        ])

    @property
    def gross(self):
        # grams carbo per plant per hour
        return np.prod([
            self.gross_CO2_umol_per_m2_s,
            self._mol_per_umol,
            self._plant_per_m2,
            self._min_step_per_sec,
            Weight.CH2O,
        ])

    @property
    def net(self):
        # grams carbo per plant per hour
        return np.prod([
            self.net_CO2_umol_per_m2_s,
            self._mol_per_umol,
            self._plant_per_m2,
            self._min_step_per_sec,
            Weight.CH2O,
        ])

    @property
    def transpiration(self):
        # Units of Transpiration from sunlit->ET are mol m-2 (leaf area) s-1
        # Calculation of transpiration from ET involves the conversion to gr per plant per hour
        #FIXME _min_step_per_sec used instead of fixed 3600 = 60 * 60
        return np.prod([
            self.transpiration_H2O_mol_per_m2_s,
            self._plant_per_m2,
            self._min_step_per_sec,
            Weight.H2O,
        ])

    @property
    def temperature(self):
        return self._weighted(self.temperature_array)

    @property
    def vapor_pressure_deficit(self):
        #HACK only use sunlit leaves?
        return np.fmax(0, self.sunlit.VPD)

    @property
    def conductance(self):
        #HACK ensure 0 when one if LAI is 0, i.e., night
        if (self.leaf_area_index_array == 0).any():
            return 0
        try:
            # average stomatal conductance Yang
            c = self._weighted(self.conductance_array) / self.p.area.leaf_area_index
            return np.fmax(0, c)
        except ZeroDivisionError:
            return 0
