from .base import Stage
from ..tracker import LeafInductionRate, BetaFunc
from .recorder import GstRecorder

import numpy as np

class TasselInitiation(Stage):
    #TODO manage juvenile/adult leaves in more general place (i.e. LeafInitiation/Appearance or Manager?)
    def setup(self, juvenile_leaves=15):
        self._juvenile_leaves = juvenile_leaves
        self._appeared_leaves_on_finish = 0
        self._temperature_recorder = GstRecorder(self.pheno)

    def tracker(self):
        return LeafInductionRate(
            pheno=self.pheno, #FIXME to access weather.day_length
            juvenile_leaves=self.juvenile_leaves,
        )

    def update(self, T):
        if self._tracker.empty():
            T = self.pheno.gst_recorder.rate
        else:
            #TODO implement on_first_update() interface?
            #HACK use mean temperature tracker for induction period
            #TODO use chained methods?
            self._temperature_recorder.update(T)
            T = self._temperature_recorder.rate
        super().update(T)

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


class PostTasselInitiation(Stage):
    def setup(self, R_max_LTAR=0.53):
        self.R_max = R_max_LTAR

    def tracker(self):
        return BetaFunc(R_max=self.R_max)

    def ready(self):
        return self.pheno.tassel_initiation.over()

    @property
    def phyllochrons_between_tassel_initiation_and_silking(self):
        return self.pheno.leaves_initiated - self.pheno.tassel_initiation._appeared_leaves_on_finish

    @property
    def rate(self):
        if self.ready():
            return super().rate / self.phyllochrons_between_tassel_initiation_and_silking
        else:
            return 0
