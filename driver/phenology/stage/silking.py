from .base import Stage
from .. import tracker

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
