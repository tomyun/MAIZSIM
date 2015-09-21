from .base import Stage
from ..tracker import GrowingDegreeDays

class GrainFillingInitiation(Stage):
    # where is this number '170' from? SK
    def setup(self, GDD_grain=170):
        self.GDD_grain = GDD_grain

    def tracker(self):
        #TODO GTI was found more accurate for grain filling stage, See Thijs phenolog paper (2014)
        return GrowingDegreeDays()

    def ready(self):
        return self.pheno.silk.over()

    def over(self):
        return self.rate >= self.GDD_grain

    def finish(self):
         GDD_sum = self.pheno.gdd_recorder.rate
         T_grow = self.pheno.gst_recorder.rate
         print("* Grain filling begins: GDDsum = {}, Growing season T = {}".format(GDD_sum, T_grow))
