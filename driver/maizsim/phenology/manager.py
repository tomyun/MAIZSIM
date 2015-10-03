from . import stage

class Phenology:
    def __init__(self, plant):
        self.plant = plant
        #HACK assumes daily timestep
        self.timestep = plant.initials.timestep / (24 * 60) # converting minute to day decimal, 1 = a day
        self.setup()

    def setup(self):
        # mean growing season temperature since germination, SK 1-19-12
        self.gst_recorder = gstr = stage.GstRecorder(self)
        self.gdd_recorder = gddr = stage.GddRecorder(self)
        self.gti_recorder = gtir = stage.GtiRecorder(self)

        self.germination = g = stage.Germination(self)
        self.emergence = e = stage.Emergence(self)
        self.leaf_initiation = li = stage.LeafInitiation(self)
        self.leaf_appearance = la = stage.LeafAppearance(self)
        self.tassel_initiation = ti = stage.TasselInitiation(self)
        self.silk = s = stage.Silk(self)
        self.grain_filling_initiation = gfi = stage.GrainFillingInitiation(self)
        self.mature = m1 = stage.Mature(self)
        self.maturity = m2 = stage.Maturity(self)
        self.death = d = stage.Death(self)

        #TODO remove PtiTracker; can be handled by LeafAppearance and TasselInitiation
        self.pti_tracker = ptit = stage.PtiTracker(self)

        self.stages = [
            gstr, gddr, gtir,
            g, e, li, la, ti, s, gfi, m1, m2, d,
            ptit,
        ]

    def __getitem__(self, index):
        return self.stages[index]

    def _queue(self):
        return [s for s in self.stages if s.ready() and not s.over()]

    def update(self):
        queue = self._queue()
        T = self.temperature
        [s.update(T) for s in queue]

        #FIXME remove finish() for simplicity
        [s.finish() for s in queue if s.over()]

    #TODO some methods for event records? or save them in Stage objects?
    #def record(self, ...):
    #    pass

    ############
    # Accessor #
    ############

    #HACK used to be leaves_total, but renamed to avoid confusion
    #HACK same as leaves_generic
    @property
    def leaves_potential(self):
        return self.plant.variety.generic_leaf_number

    @property
    def leaves_generic(self):
        return self.plant.variety.generic_leaf_number

    @property
    def leaves_initiated(self):
        return self.leaf_initiation.leaves

    @property
    def leaves_appeared(self):
        return self.leaf_appearance.leaves

    @property
    def temperature(self):
        if self.leaves_appeared < 9:
            T = self.plant.soil.T_soil
        else:
            T = self.plant.weather.T_air
        #FIXME T_cur doesn't go below zero, but is it fair assumption?
        return T if T > 0 else 0

    @property
    def growing_temperature(self):
        return self.gst_recorder.rate

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
    def tassel_initiated(self):
        return self.tassel_initiation.over()

    @property
    def vegetative_growing(self):
        #HACK there would be a missing period right after tassel initiation and before silking, waiting for full leaf appearance
        #return self.germinated and not self.tassel_initiated
        return not self.silking

    @property
    def silking(self):
        return self.silk.ing()

    @property
    def silked(self):
        return self.silk.over()

    @property
    def grain_filling(self):
        return self.grain_filling_initiation.over()

    @property
    def matured(self):
        return self.mature.over()

    @property
    def dead(self):
        return self.death.over()

    # GDDsum
    @property
    def gdd_after_emergence(self):
        if self.emergence.over():
            #HACK tracker is reset when emergence is over
            return self.gst_recorder.rate
        else:
            return 0

    #TODO serve them directly from TasselInitiation?
    @property
    def leaves_appeared_since_tassel_initiation(self):
        #FIXME no need to use a separate tracker
        #return self.pti_tracker.rate
        if self.tassel_initiation.over():
            return self.leaves_appeared - self.tassel_initiation.appeared_leaves_on_finish
        else:
            return 0

    @property
    def leaves_to_appear_since_tassel_initiation(self):
        if self.tassel_initiation.over():
            return self.leaves_initiated - self.tassel_initiation.appeared_leaves_on_finish
        else:
            return 0

    @property
    def leaf_appearance_fraction_since_tassel_initiation(self):
        if self.tassel_initiation.over():
            return self.leaves_appeared_since_tassel_initiation / self.leaves_to_appear_since_tassel_initiation
        else:
            return 0

    @property
    def current_stage(self):
        if self.matured:
            return "Matured"
        elif self.grain_filling:
            return "grainFill"
        elif self.silked:
            return "Silked"
        #FIXME no flowered (anthesis)?
        #elif self.flowered:
            #return "Flowered"
        elif self.tassel_initiated:
            return "Tasselinit"
        elif self.emerged:
            return "Emerged"
        elif self.dead:
            return "Inactive"
        else:
            return "none"
