from organ import Organ

class Stem(Organ):
    def __init__(self, nodal_unit):
        self.nodal_unit = nodal_unit

    def setup(self):
        self.rank = self.nodal_unit.rank
        self.length = 0
        self.diameter = 0

    def update(self):
        super().update(self.pheno.temperature)
