from .base import Stage
from .. import tracker

class LeafInitiation(Stage):
    #FIXME use correct args
    def setup(self, initial_leaves=5, R_max_LIR=0.0978):
        self.initial_leaves = initial_leaves
        self.leaves = initial_leaves
        self.R_max = R_max_LIR

    def tracker(self):
        t = tracker.BetaFunc(self.R_max)
        t.set_initial_value(self.initial_leaves)
        return t

    def post_update(self):
        self.leaves = self.initial_leaves + int(self.rate)

    def ready(self):
        return self.pheno.germination.over()

    def over(self):
        return self.pheno.tassel_initiation.over()
