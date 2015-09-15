#from ..tracker import BetaFunc, GrowingDegreeDays, VegetativeGeneralThermalIndex, ReproductiveGeneralThermalIndex, LeafInductionRate
from .. import tracker

import numpy as np

class Stage(object):
    def __init__(self, pheno):
        self.pheno = pheno
        self.setup()
        self._tracker = self.tracker()

    def setup(self):
        pass

    #TODO use @ decorator instead?
    def tracker(self):
        raise NotImplementedError("Need to create a list of Tracker objects.")

    #TODO prevent duplicate updates
    def update(self, T):
        dt = self.pheno.timestep # per day
        self._tracker.update(T, dt)

    def post_update(self):
        pass

    @property
    def rate(self):
        return self._tracker.rate

    #TODO make ready, over, ing propery?

    def ready(self):
        return False

    def over(self):
        return False

    def ing(self):
        return self.ready() and not self.over()

    #HACK better way to handle callback?
    def finish(self):
        pass
