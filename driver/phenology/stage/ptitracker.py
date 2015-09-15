from .base import Stage
from ..tracker import BetaFunc

#FIXME better naming... isn't it a duplicate of Silking?
# to be used for C partitoining time scaling, see Plant.cpp
class PtiTracker(Stage):
    #FIXME use correct args
    def setup(self, R_max_LTAR=0.53):
        self.R_max = R_max_LTAR

    def tracker(self):
        return BetaFunc(R_max=self.R_max)

    def ready(self):
        return self.pheno.tassel_initiation.over()
