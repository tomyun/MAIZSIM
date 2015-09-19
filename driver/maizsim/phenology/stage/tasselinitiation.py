from .base import Stage
from ..tracker import LeafInductionRate

import numpy as np

class TasselInitiation(Stage):
    #FIXME use correct args
    def setup(self, juvenile_leaves=15, day_length=None):
        self.juvenile_leaves = juvenile_leaves
        self.day_length = day_length
        self._appeared_leaves = None

    def tracker(self):
        return LeafInductionRate(
            pheno=self.pheno, #FIXME to access weather.day_length
            gst_tracker=self.pheno.gst_tracker,
            juvenile_leaves=self.juvenile_leaves,
        )

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
