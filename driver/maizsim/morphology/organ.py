from .weight import Weight
from ..phenology.tracker import GrowingDegreeDays

class Organ:
    def __init__(self, plant):
        self.p = plant
        self._tracker = GrowingDegreeDays(T_base=8.0, T_opt=None, T_max=43.3).use_timestep(plant.pheno.timestep)

        # organ temperature, C
        self._temperature = 25.0

        # glucose, MW = 180.18 / 6 = 30.03 g
        self._carbohydrate = 0

        # nitrogen content, mg
        self._nitrogen = 0

        self.setup()

    def setup(self):
        pass

    # chronological age of an organ, days
    @property
    def age(self):
        return self._tracker.period

    # physiological age accouting for temperature effect (in reference to endGrowth and lifeSpan, days)
    @property
    def physiological_age(self):
        return self._tracker.rate

    # biomass, g
    @property
    def mass(self):
        #FIXME isn't it just the amount of carbohydrate?
        #return self._carbohydrate / Weight.CH2O * Weight.C / Weight.C_to_CH2O_ratio
        return self._carbohydrate

    #TODO remove set_mass() and directly access carbohydrate
    def set_mass(self, mass):
        #self._carbohydrate = mass * Weight.C_to_CH2O_ratio / Weight.C * Weight.CH2O
        self._carbohydrate = mass

    # physiological days to reach the end of growth (both cell division and expansion) at optimal temperature, days
    @property
    def growth_duration(self):
        return 10

    # life expectancy of an organ in days at optimal temperature (fastest growing temp), days
    #FIXME not used
    @property
    def longevity(self):
        return 50

    # carbon allocation to roots or leaves for time increment
    #FIXME not used
    @property
    def potential_carbohydrate_increment(self):
        return 0

    # carbon allocation to roots or leaves for time increment  gr C for roots, gr carbo dt-1
    #FIXME not used
    @property
    def actual_carbohydrate_increment(self):
        return 0

    def _update(self, T):
        self._temperature = T
        self._tracker.update(T)

    def update(self):
        self._update(self.p.pheno.temperature)

    def import_carbohydrate(self, amount):
        self._carbohydrate += amount

    def import_nitrogen(self, amount):
        self._nitrogen += amount

    def respire(self):
        # this needs to be worked on
        # currently not used at all
        Ka = 0.1 # growth respiration
        Rm = 0.02 # maintenance respiration
        #self._carbohydrate -= (Ka + Rm) * self._carbohydrate
        self._carbohydrate *= 1 - (Ka + Rm)
