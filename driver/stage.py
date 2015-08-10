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
        raise NotImplementedError("Need to create a list of ThermalFunc objects.")

    #TODO prevent duplicate updates
    def update(self, T):
        dt = self.pheno.timestep # per day
        self._tracker.update(T, dt)

    def post_update(self):
        pass

    @property
    def rate(self):
        return self._tracker.rate

    def ready(self):
        return False

    def over(self):
        return False

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
    def setup(self, R_max_LIR):
        self.R_max = R_max_LIR
        self.leaves = 0

    def tracker(self):
        return tracker.BetaFunc(self.R_max)

    def post_update(self):
        self.leaves = int(self.rate)

    def ready(self):
        return self.pheno.germination.over()

    def over(self):
        return self.pheno.tassel_initiation.over()


class TasselInitiation(Stage):
    def setup(self, juvenile_leaves, day_length):
        self.juvenile_leaves = juvenile_leaves
        self.day_length = day_length

    def tracker(self):
        return tracker.LeafInductionRate(self.pheno.gst_tracker, self.juvenile_leaves, self.day_length)

    @property
    def initiated_leaves(self):
        return self.pheno.leaf_initiation.leaves

    @property
    def added_leaves(self):
        return self.initiated_leaves - self.juvenile_leaves

    @property
    def leaves_to_induce(self):
        return int(self.rate)

    def ready(self):
        return self.added_leaves >= 0

    def over(self):
        return self.added_leaves >= self.leaves_to_induce

    def finish(self):
        #TODO clean up leaf count variables
        self.current_leaf = self.youngest_leaf = self.total_leaves = self.initiated_leaves
        #self.leaves_at_tassel_initiation = leaves_appeared?

        GDD_sum = self.pheno.gdd_tracker.rate
        T_grow = self.pheno.gst_tracker.rate
        print("* Tassel initiation: GDDsum = {}, Growing season T = {}".format(GDD_sum, T_grow))


#FIXME better naming
# to be used for C partitoining time scaling, see Plant.cpp
class PhyllochronsFromTI(Stage):
    def setup(self, R_max_LTAR):
        self.R_max = R_max_LTAR

    def tracker(self):
        return tracker.BetaFunc(self.R_max)

    def ready(self):
        return self.pheno.tassel_initiation.over()


class LeafAppearance(Stage):
    def setup(self, R_max_LTAR):
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
    def setup(self, R_max_LTAR, phyllochrons):
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


class GrainFilling(Stage):
    # where is this number '170' from? SK
    def setup(self, GDD_grain=170):
        self.GDD_grain = GDD_grain

    def tracker(self):
        #TODO GTI was found more accurate for grain filling stage, See Thijs phenolog paper (2014)
        return tracker.GrowingDegreeDays()

    def ready(self):
        return self.pheno.germination.over()

    def over(self):
        return self.rate >= self.GDD_grain

    def finish(self):
         GDD_sum = self.pheno.gdd_tracker.rate
         T_grow = self.pheno.gst_tracker.rate
         print("* Grain filling begins: GDDsum = {}, Growing season T = {}".format(GDD_sum, T_grow))


class GddTracker(Stage):
    def create_trackers(self):
        return [GrowingDegreeDays()]


class GtiTracker(Stage):
    def create_trackers(self):
        return [ReproductiveGeneralThermalIndex()]
