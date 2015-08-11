from tracker import GrowingDegreeDays

class Weight(object):
    CO2 = 44.0098
    C = 12.011
    CH2O = 30.03


class Organ(object):
    def __init__(self):
        self._tracker = GrowingDegreeDays(timestep=1/24/60, T_base=8.0, T_opt=None, T_max=43.3)

        # organ temperature, C
        self.temperature = 25.0

        # glucose, MW = 180.18 / 6 = 30.03 g
        self.carbohydrate = 0

        # nitrogen content, mg
        self.nitrogen = 0

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
        #return self.carbohydrate
        C_to_CH2O_ratio = Weight.C / Weight.CH2O # 0.40
        return self.carbohydrate / Weight.CH2O * Weight.C / C_to_CH2O_ratio

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

    def update(self, T):
        self.temperature = T
        self._tracker.update(T)

    def import_carbohydrate(self, amount):
        self.carbohydrate += amount

    def import_nitrogen(self, amount):
        self.nitrogen += amount

    def respire(self):
        # this needs to be worked on
        # currently not used at all
        Ka = 0.1 # growth respiration
        Rm = 0.02 # maintenance respiration
        #self.carbohydrate -= (Ka + Rm) * self.carbohydrate
        self.carbohydrate *= 1 - (Ka + Rm)
