from . import stage

class Phenology(object):
    def __init__(self, plant):
        self.plant = plant
        #TODO make use of InitInfo object
        self.timestep = plant.info.timestep
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
        self.mature = m1 = stage.Mature(self)
        self.maturity = m2 = stage.Maturity(self)
        self.death = d = stage.Death(self)

        #TODO remove PtiTracker; can be handled by LeafAppearance and TasselInitiation
        self.pti_tracker = ptit = stage.PtiTracker(self)

        self.stages = [
            gstt, gddt, gtit,
            g, e, li, la, ti, s, gfi, m1, m2, d,
            ptit,
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

    #TODO some methods for event records? or save them in Stage objects?
    #def record(self, ...):
    #    pass

    ############
    # Accessor #
    ############

    #FIXME same as leaves_initiated
    @property
    def leaves_total(self):
        return self.leaf_initiation.leaves

    @property
    def leaves_generic(self):
        #TODO from TInitInfo
        return 15

    @property
    def leaves_initiated(self):
        return self.leaf_initiation.leaves

    @property
    def leaves_appeared(self):
        return self.leaf_appearance.leaves

    #TODO is it relevant here?
    @property
    def temperature(self):
        #TODO modify update() to get Atmos object...
        return T

    @property
    def growing_temperature(self):
        return self.gst_tracker.rate

    @property
    def optimal_temperature(self):
        #TODO parmaterize?
        return 32.1

    @property
    def germinating(self):
        return self.germination.ing()

    @property
    def germinated(self):
        return self.germination.over()

    @property
    def emerging(self):
        return self.emergence.ing()

    @property
    def emerged(self):
        return self.emergence.over()

    @property
    def vegetative_growing(self):
        return self.germination.over() and not self.tassel_initiation.over()

    @property
    def silking(self):
        return self.silking.ing()

    @property
    def grain_filling(self):
        return self.grain_filling_initiation.over()

    @property
    def dead(self):
        return self.death.over()

    @property
    def gdd_after_emergence(self):
        if self.emergence.over():
            #HACK tracker is reset when emergence is over
            return self.gst_tracker.rate
        else:
            return 0

    @property
    def leaves_appeared_since_tassel_initiation(self):
        #FIXME no need to use a separate tracker
        #return self.pti_tracker.rate
        return self.leaves_appeared - self.tassel_initiation.appeared_leaves

    @property
    def leaves_to_appear_since_tassel_initiation(self):
        return self.tassel_initiation.leaves_to_appear

    @property
    def leaf_appearance_fraction_since_tassel_initiation(self):
        return self.leaves_appeared_since_tassel_initiation / self.leaves_to_appear_since_tassel_initiation
