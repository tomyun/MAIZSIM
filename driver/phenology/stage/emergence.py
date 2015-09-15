from .base import Stage
from .. import tracker

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
