from .base import Stage
from ..tracker import BetaFunc

class LeafInitiation(Stage):
    #FIXME use correct args
    def setup(self, initial_leaves=5, R_max_LIR=0.978):
        self.initial_leaves = initial_leaves
        self.R_max = R_max_LIR

    def tracker(self):
        return BetaFunc(R_max=self.R_max)

    def ready(self):
        return self.pheno.germination.over()

    def over(self):
        return self.pheno.tassel_initiation.over()

    @property
    def leaves(self):
        return self.initial_leaves + int(self.rate)
