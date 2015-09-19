from .trait import Trait
from ..morphology.weight import Weight

import numpy as np

#TODO rename to CarbonAssimilation or so? could be consistently named as CarbonPartition, CarbonAllocation...
class Photosynthesis(Trait):
    def setup(self):
        self.sunlit_leaf = GasExchange('Sunlit')
        self.shaded_leaf = GasExchange('Shaded')

    @property
    def sunlit_leaf_area_index(self):
        return self.p.lightenv.sunlitLAI()

    @property
    def shaded_leaf_area_index(self):
        return self.p.lightenv.shadedLAI()

    @property
    def leaf_area_index_array(self):
        return np.array([
            self.sunlit_leaf_area_index,
            self.shaded_leaf_area_index,
        ])

    def _weighted(self, array):
        return self.leaf_area_index_array.dot(array)

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
        return 1 / initInfo.plant_density

    @property
    def _min_step_per_sec(self):
        return 60 * initInfo.time_step

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
        #TODO is this condition necessary?
        #if sunlit_LAI >= 0 and shaded_LAI >= 0 and LAI >= 0:
        try:
            # average stomatal conductance Yang
            c = self._weighted(self.conductance_array) / self.p.area.leaf_area_index
            return np.fmax(0, c)
        except ZeroDivisionError:
            return 0
