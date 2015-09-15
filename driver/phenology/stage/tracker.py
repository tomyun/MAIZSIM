from .base import Stage
from .. import tracker

# Non-growth related Stage classes for tracking thermal units over entire growth period
class TrackerStage(Stage):
    def ready(self):
        return True

    def reset(self):
        self._tracker.reset()


class GstTracker(TrackerStage):
    def tracker(self):
        return tracker.Tracker()


class GddTracker(TrackerStage):
    def tracker(self):
        return tracker.GrowingDegreeDays()


class GtiTracker(TrackerStage):
    def tracker(self):
        return tracker.ReproductiveGeneralThermalIndex()
