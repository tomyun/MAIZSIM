#from .tracker import BetaFunc, GrowingDegreeDays, VegetativeGeneralThermalIndex, ReproductiveGeneralThermalIndex, LeafInductionRate
import tracker

import numpy as np

class Stage(object):
    def __init__(self, pheno):
        self.pheno = pheno
        self.setup()
        self._tracker = self.tracker()

    def setup(self):
        pass

    #TODO use @ decorator instead?
    def tracker(self):
        raise NotImplementedError("Need to create a list of Tracker objects.")

    #TODO prevent duplicate updates
    def update(self, T):
        dt = self.pheno.timestep # per day
        self._tracker.update(T, dt)

    def post_update(self):
        pass

    @property
    def rate(self):
        return self._tracker.rate

    #TODO make ready, over, ing propery?

    def ready(self):
        return False

    def over(self):
        return False

    def ing(self):
        return self.ready() and not self.over()

    #HACK better way to handle callback?
    def finish(self):
        pass


class Germination(Stage):
    def setup(self, R_max=0.45):
        self.R_max = R_max

    def tracker(self):
        return tracker.BetaFunc(self.R_max)

    def ready(self):
        #TODO implement germination rate model of temperature
        # for now assume it germinates immidiately after sowing
        return True

    def over(self):
        return self.rate >= 0.5

    def finish(self):
        GDD_sum = self.pheno.gdd_tracker.rate
        dt = self.pheno.timestep * 24 * 60 # per min
        print("* Germinated: GDDsum = {}, time step (min) = {}".format(GDD_sum, dt))


class Emergence(Stage):
    def setup(self, R_max=0.2388):
        self.R_max = R_max

    def tracker(self):
        return tracker.BetaFunc(self.R_max)

    def ready(self):
        return self.pheno.germination.over()

    def over(self):
        return self.rate >= 1.0

    def finish(self):
        GDD_sum = self.pheno.gdd_tracker.rate
        T_grow = self.pheno.gst_tracker.rate
        print("* Emergence: GDDsum = {}, Growing season T = {}".format(GDD_sum, T_grow))

        #HACK reset GDD tracker after emergence
        self.emerge_GDD = GDD_sum
        self.pheno.gdd_tracker.reset()


class LeafInitiation(Stage):
    #FIXME use correct args
    def setup(self, initial_leaves=5, R_max_LIR=0.0978):
        self.initial_leaves = initial_leaves
        self.leaves = initial_leaves
        self.R_max = R_max_LIR

    def tracker(self):
        t = tracker.BetaFunc(self.R_max)
        t.set_initial_value(self.initial_leaves)
        return t

    def post_update(self):
        self.leaves = self.initial_leaves + int(self.rate)

    def ready(self):
        return self.pheno.germination.over()

    def over(self):
        return self.pheno.tassel_initiation.over()


class TasselInitiation(Stage):
    #FIXME use correct args
    def setup(self, juvenile_leaves=15, day_length=None):
        self.juvenile_leaves = juvenile_leaves
        self.day_length = day_length
        self._appeared_leaves = None

    def tracker(self):
        return tracker.LeafInductionRate(self.pheno.gst_tracker, self.juvenile_leaves, self.day_length)

    @property
    def initiated_leaves(self):
        return self.pheno.leaf_initiation.leaves

    @property
    def added_leaves(self):
        return np.fmin(0, self.initiated_leaves - self.juvenile_leaves)

    @property
    def leaves_to_induce(self):
        return self.rate

    @property
    def appeared_leaves(self):
        return self._appeared_leaves

    @property
    def leaves_to_appear(self):
        return self.initiated_leaves - self.appeared_leaves

    def ready(self):
        return self.added_leaves >= 0

    def over(self):
        return self.added_leaves >= self.leaves_to_induce

    def finish(self):
        #TODO clean up leaf count variables
        self.current_leaf = self.youngest_leaf = self.total_leaves = self.initiated_leaves
        #HACK save the appeared leaves when tassel initiation is done
        self._appeared_leaves = self.pheno.leaf_appearance.leaves

        GDD_sum = self.pheno.gdd_tracker.rate
        T_grow = self.pheno.gst_tracker.rate
        print("* Tassel initiation: GDDsum = {}, Growing season T = {}".format(GDD_sum, T_grow))


#FIXME better naming... isn't it a duplicate of Silking?
# to be used for C partitoining time scaling, see Plant.cpp
class PtiTracker(Stage):
    #FIXME use correct args
    def setup(self, R_max_LTAR=0.53):
        self.R_max = R_max_LTAR

    def tracker(self):
        return tracker.BetaFunc(self.R_max)

    def ready(self):
        return self.pheno.tassel_initiation.over()


class LeafAppearance(Stage):
    #FIXME use correct args
    def setup(self, R_max_LTAR=0.53):
        self.R_max = R_max_LTAR
        self.leaves = 0

    def tracker(self):
        return tracker.BetaFunc(self.R_max)

    def post_update(self):
        self.leaves = int(self.rate)

    def ready(self):
        initiated_leaves = self.pheno.leaf_initiation.leaves
        return self.leaves < initiated_leaves

    def over(self):
        initiated_leaves = self.pheno.leaf_initiation.leaves
        #HACK ensure leaves are initiated
        return self.leaves >= initiated_leaves > 0


class Silking(Stage):
    #FIXME use correct args
    #TODO check the correct phyllochrons: is it 8 or 3?
    def setup(self, R_max_LTAR=0.53, phyllochrons=8):
        self.R_max = R_max_LTAR
        self.phyllochrons = phyllochrons

    def tracker(self):
        # Assume 75% Silking occurs at total tip appeared + 3 phyllochrons
        return tracker.BetaFunc(self.R_max)

    def ready(self):
        return self.pheno.tassel_initiation.over() and self.pheno.leaf_appearance.over()

    def over(self):
        # anthesis rate
        return self.rate >= self.phyllochrons

    def finish(self):
         GDD_sum = self.pheno.gdd_tracker.rate
         T_grow = self.pheno.gst_tracker.rate
         print("* Silking: GDDsum = {}, Growing season T = {}".format(GDD_sum, T_grow))


class GrainFillingInitiation(Stage):
    # where is this number '170' from? SK
    def setup(self, GDD_grain=170):
        self.GDD_grain = GDD_grain

    def tracker(self):
        #TODO GTI was found more accurate for grain filling stage, See Thijs phenolog paper (2014)
        return tracker.GrowingDegreeDays()

    def ready(self):
        return self.pheno.silking.over()

    def over(self):
        return self.rate >= self.GDD_grain

    def finish(self):
         GDD_sum = self.pheno.gdd_tracker.rate
         T_grow = self.pheno.gst_tracker.rate
         print("* Grain filling begins: GDDsum = {}, Growing season T = {}".format(GDD_sum, T_grow))


class Mature(Stage):
    def setup(self, GDD_rating=1331):
        self.GDD_rating = GDD_rating

    def tracker(self):
        return tracker.GrowingDegreeDays()

    def ready(self):
        return self.pheno.emergence.over()

    def over(self):
        return self.rate >= self.GDD_rating

    def finish(self):
         GDD_sum = self.pheno.gdd_tracker.rate
         T_grow = self.pheno.gst_tracker.rate
         print("* Matured: rate = {}, GDDsum = {}, Growing season T = {}".format(self.rate, GDD_sum, T_grow))


#FIXME quite confusing names: Mature vs. Maturity
class Maturity(Stage):
    def tracker(self):
        #HACK no use
        return tracker.Tracker()

    def ready(self):
        return True

    def over(self):
        # when less than 5% of total leaf area is green, physiological maturity is reached. SK
        # see http://www.agry.purdue.edu/ext/corn/news/timeless/TopLeafDeath.html
        area = self.pheno.plant.area
        return area.green_leaf <= 0.05 * area.leaf

    def finish(self):
        GDD_sum = self.pheno.gdd_tracker.rate
        active_leaf_percent = self.pheno.plant.area.active_leaf_ratio * 100
        print("* Physiological maturity: GDDsum = {}, Growing season T = {}, green leaf: {}%".format(GDD_sum, T_grow, active_leaf_percent))


class Death(Stage):
    def tracker(self):
        #HACK no use
        return tracker.Tracker()

    def ready(self):
        return True

    def over(self):
        #HACK access pheno and plant from stage...
        return self.pheno.plant.count.total_dropped_leaves >= self.pheno.leaves_total

    def finish(self):
        #FIXME record event?
        GDD_sum = self.pheno.gdd_tracker.rate
        print("* Death: GDDsum = {}, Growing season T = {}".format(GDD_sum, T_grow))


# Non-growth related Stage classes for tracking thermal units over entire growth period
class TrackerStage(Stage):
    def ready(self):
        return True

    def reset(self):
        self._tracker.reset()


class GstTracker(TrackerStage):
    def tracker(self):
        return tracker.Tracker()


class GddTracker(TrackerStage):
    def tracker(self):
        return tracker.GrowingDegreeDays()


class GtiTracker(TrackerStage):
    def tracker(self):
        return tracker.ReproductiveGeneralThermalIndex()
