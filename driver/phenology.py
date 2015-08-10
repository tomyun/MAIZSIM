import numpy as np

import stage

class Phenology(object):
    def __init__(self, timestep):
        self.timestep = timestep
        self.setup()

    def setup(self):
        # mean growing season temperature since germination, SK 1-19-12
        self.gst_tracker = gstt = stage.GstTracker(self)
        self.gdd_tracker = gddt = stage.GddTracker(self)
        self.gti_tracker = gtit = stage.GtiTracker(self)

        self.germination = g = stage.Germination(self)
        self.emergence = e = stage.Emergence(self)
        self.leaf_initiation = li = stage.LeafInitiation(self)
        self.leaf_appearance = la = stage.LeafAppearance(self)
        self.tassel_initiation = ti = stage.TasselInitiation(self)
        self.silking = s = stage.Silking(self)
        self.grain_filling_initiation = gfi = stage.GrainFillingInitiation(self)
        self.mature = m = stage.Mature(self)
        #self.maturity = m = Maturity(self)

        self.phyllochrons_from_ti = pti = stage.PhyllochronsFromTI(self)

        self.stages = [
            gstt, gddt, gtit,
            g, e, li, la, ti, s, gfi, m,
            pti,
        ]

    def __getitem__(self, index):
        return self.stages[index]

    def _queue(self):
        return [s for s in self.stages if s.ready() and not s.over()]

    def update(self, T):
        queue = self._queue()

        [s.update(T) for s in queue]
        [s.post_update() for s in queue]

        #FIXME remove finish() for simplicity
        [s.finish() for s in queue if s.over()]
