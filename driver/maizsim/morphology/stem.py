from .organ import Organ

class Stem(Organ):
    def __init__(self, nodal_unit):
        super().__init__(nodal_unit.plant)
        self.nodal_unit = nodal_unit

    def setup(self):
        self.length = 0
        self.diameter = 0

    @property
    def rank(self):
        return self.nodal_unit.rank

    def update(self):
        super().update()
