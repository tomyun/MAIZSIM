from .base import Stage
from ..tracker import BetaFunc

class Germination(Stage):
    def setup(self, R_max=0.45):
        self.R_max = R_max

    def tracker(self):
        return BetaFunc(R_max=self.R_max)

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
