from .base import Stage
from ..tracker import LeafInductionRate
from .recorder import GstRecorder

import numpy as np

class TasselInitiation(Stage):
    #TODO manage juvenile/adult leaves in more general place (i.e. LeafInitiation/Appearance or Manager?)
    def setup(self, juvenile_leaves=15, day_length=None):
        self._juvenile_leaves = juvenile_leaves
        self.day_length = day_length
        self._appeared_leaves_on_finish = 0

    def tracker(self):
        return LeafInductionRate(
            pheno=self.pheno, #FIXME to access weather.day_length
            gst_recorder=self.pheno.gst_recorder,
            temperature_recorder=GstRecorder(self.pheno),
            juvenile_leaves=self.juvenile_leaves,
        )

    @property
    def juvenile_leaves(self):
        return self._juvenile_leaves

    @property
    def _adult_leaves(self):
        return self.pheno.leaves_initiated - self.juvenile_leaves

    @property
    def adult_leaves(self):
        return np.fmax(0, self._adult_leaves)

    @property
    def leaves_to_induce(self):
        return self.rate

    #FIXME is it really needed?
    @property
    def appeared_leaves_on_finish(self):
        return self._appeared_leaves_on_finish

    def ready(self):
        #HACK should be negative when not ready
        return self._adult_leaves >= 0

    def over(self):
        return self.adult_leaves >= self.leaves_to_induce

    def finish(self):
        #TODO clean up leaf count variables
        self.current_leaf = self.youngest_leaf = self.total_leaves = self.pheno.leaves_initiated
        #HACK save the appeared leaves when tassel initiation is done
        self._appeared_leaves_on_finish = self.pheno.leaves_appeared

        GDD_sum = self.pheno.gdd_recorder.rate
        T_grow = self.pheno.gst_recorder.rate
        print("* Tassel initiation: GDDsum = {}, Growing season T = {}".format(GDD_sum, T_grow))
