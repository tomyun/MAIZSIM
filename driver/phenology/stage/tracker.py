from .base import Stage
from ..tracker import Tracker, GrowingDegreeDays, ReproductiveGeneralThermalIndex

# Non-growth related Stage classes for tracking thermal units over entire growth period
class TrackerStage(Stage):
    def ready(self):
        return True

    def reset(self):
        self._tracker.reset()


class GstTracker(TrackerStage):
    def tracker(self):
        return Tracker()


class GddTracker(TrackerStage):
    def tracker(self):
        return GrowingDegreeDays()


class GtiTracker(TrackerStage):
    def tracker(self):
        return ReproductiveGeneralThermalIndex()
