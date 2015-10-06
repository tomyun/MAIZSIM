from .base import Stage
from ..tracker import TemperatureTracker, GrowingDegreeDays, ReproductiveGeneralThermalIndex

# Non-growth related Stage classes for tracking thermal units over entire growth period
class Recorder(Stage):
    def ready(self):
        return True

    def reset(self):
        self._tracker.reset()


class GstRecorder(Recorder):
    def tracker(self):
        return TemperatureTracker()

    #HACK growing season only counts after germination
    def ready(self):
        return self.pheno.germination.over()


class GddRecorder(Recorder):
    def tracker(self):
        return GrowingDegreeDays()


class GtiRecorder(Recorder):
    def tracker(self):
        return ReproductiveGeneralThermalIndex()
