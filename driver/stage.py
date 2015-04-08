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
        return BetaFunc(self.R_max)

    def ready(self):
        #TODO implement germination rate model of temperature
        # for now assume it germinates immidiately after sowing
        return True

    def over(self):
        return self.rate >= 0.5

    def finish(self):
        GDD_sum = self.pheno.gdd_tracker.rate
        dt = self.pheno.timestep * 24 * 60 # per min
        print("* Germinated: GDDsum = {}, time step (min) = {}" % (GDD_sum, dt))


class Emergence(Stage):
    def setup(self, R_max=0.2388):
        self.R_max = R_max

    def tracker(self):
        return BetaFunc(self.R_max)

    def ready(self):
        return self.pheno.germination.over()

    def over(self):
        return self.rate >= 1.0

    def finish(self):
        GDD_sum = self.pheno.gdd_tracker.rate
        T_grow = self.pheno.gst_tracker.rate
        print("* Emergence: GDDsum = {}, Growing season T = {}" % (GDD_sum, T_grow))

        #HACK reset GDD tracker after emergence
        self.emerge_GDD = GDD_sum
        self.pheno.gdd_tracker.reset()


class LeafInitiation(Stage):
    def setup(self, R_max_LIR):
        self.R_max = R_max_LIR

    def tracker(self):
        return BetaFunc(self.R_max)

    def post_update(self):
        self.leaves = int(self.rate)

    def ready(self):
        return self.pheno.germination.over()

    #FIXME should keep going until tassel initiation
    def over(self):
        return self.pheno.tassel_initiation.ready()


class TasselInitiation(Stage):
    def setup(self, leaves_to_induce, juvenile_leaves, day_length):
        self.leaves_to_induce = leaves_to_induce
        self.juvenile_leaves = juvenile_leaves
        self.day_length = day_length

    def tracker(self):
        #HACK create its own temperature tracker
        tt = tracker.Tracker()
        tt.update(self.pheno.gst_tracker.rate)
        return tracker.LeafInductionRate(tt, self.juvenile_leaves, self.day_length)

    def ready(self):
        initiated_leaves = self.pheno.leaf_initiation.leaves
        return initiated_leaves >= self.juvenile_leaves

    def over(self):
        initiated_leaves = self.pheno.leaf_initiation.leaves
        added_leaves = initiated_leaves - self.juvenile_leaves
        return added_leaves >= self.leaves_to_induce

    def finish(self):
        #TODO clean up leaf count variables
        initiated_leaves = self.pheno.leaf_initiation.leaves
        self.current_leaf = self.youngest_leaf = self.total_leaves = initiated_leaves
        #self.leaves_at_tassel_initiation = leaves_appeared?

        GDD_sum = self.pheno.gdd_tracker.rate
        T_grow = self.pheno.gst_tracker.rate
        print("* Tassel initiation: GDDsum = {}, Growing season T = {}" % (GDD_sum, T_grow))


#FIXME better naming
class PhyllochronsFromTI(Stage):
    pass


class LeafAppearance(Stage):
    pass


#FIXME probably indicates Silking
class AfterTasselInitiation(Stage):
    pass


class Silking(Stage):
    pass


class GrainFilling(Stage):
    pass


class GddTracker(Stage):
    def create_trackers(self):
        return [GrowingDegreeDays()]


class GtiTracker(Stage):
    def create_trackers(self):
        return [ReproductiveGeneralThermalIndex()]
