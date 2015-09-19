import numpy as np

class Stage(object):
    def __init__(self, pheno):
        self.pheno = pheno
        self.setup()
        self._tracker = self.tracker().use_timestep(pheno.timestep)

    def setup(self):
        pass

    def tracker(self):
        raise NotImplementedError("Need to create a list of Tracker objects.")

    #TODO prevent duplicate updates
    def update(self, T):
        self._tracker.update(T)

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
