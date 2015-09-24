from .base import Accumulator

import numpy as np

#TODO it is not really about phenology, better name? place?
class WaterStress(Accumulator):
    def setup(self, scale=1):
        # scale for reduction in leaf lifespan and aging rates
        self.scale = scale

    def calc(self, effect):
        # This assumes 0.25mg/m2 minimum N required, and below this the value is 0.0.
        # threshold predawn leaf water potential (in bars) below which water stress triggers senescence, needs to be substantiated with lit or exp evidence, SK
        # This is the water potential at which considerable reduction in leaf growth takes place in corn, sunflower, and soybean in Boyear (1970)
        return self.scale * (1 - effect)
