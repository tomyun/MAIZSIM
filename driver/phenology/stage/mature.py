from .base import Stage
from ..tracker import GrowingDegreeDays

class Mature(Stage):
    def setup(self, GDD_rating=1331):
        self.GDD_rating = GDD_rating

    def tracker(self):
        return GrowingDegreeDays()

    def ready(self):
        return self.pheno.emergence.over()

    def over(self):
        return self.rate >= self.GDD_rating

    def finish(self):
         GDD_sum = self.pheno.gdd_tracker.rate
         T_grow = self.pheno.gst_tracker.rate
         print("* Matured: rate = {}, GDDsum = {}, Growing season T = {}".format(self.rate, GDD_sum, T_grow))
