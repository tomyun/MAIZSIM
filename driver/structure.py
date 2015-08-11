class Stem(object):
    pass


from enum import Enum

class State(Enum):
    initiated = 1
    appeared = 2
    growing = 3
    prolific = 4
    aging = 5
    terminated = 6

class NodalUnit(object):
    def __init__(self, pheno, rank, mass):
        self.pheno = pheno
        self.rank = rank
        self.leaf = Leaf(self)
        self.stem = Stem(self)

    # @property
    # def initiated(self):
    #     return True

    @property
    def mass(self):
        return self.leaf.mass + self.stem.mass

    def update(self, predawn_lwp):
        self.leaf.update()
        self.stem.update()
