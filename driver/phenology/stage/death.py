from .base import Stage
from .. import tracker

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
