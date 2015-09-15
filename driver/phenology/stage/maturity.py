from .base import Stage
from .. import tracker

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
