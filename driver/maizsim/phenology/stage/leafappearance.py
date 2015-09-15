from .base import Stage
from ..tracker import BetaFunc

class LeafAppearance(Stage):
    #FIXME use correct args
    def setup(self, R_max_LTAR=0.53):
        self.R_max = R_max_LTAR
        self.leaves = 0

    def tracker(self):
        return BetaFunc(R_max=self.R_max)

    def post_update(self):
        self.leaves = int(self.rate)

    def ready(self):
        initiated_leaves = self.pheno.leaf_initiation.leaves
        return self.leaves < initiated_leaves

    def over(self):
        initiated_leaves = self.pheno.leaf_initiation.leaves
        #HACK ensure leaves are initiated
        return self.leaves >= initiated_leaves > 0