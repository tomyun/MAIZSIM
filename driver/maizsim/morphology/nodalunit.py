from .leaf import Leaf
from .stem import Stem

from enum import Enum

class State(Enum):
    initiated = 1
    appeared = 2
    growing = 3
    prolific = 4
    aging = 5
    terminated = 6

class NodalUnit:
    def __init__(self, plant, rank):
        self.plant = plant
        self.rank = rank
        self.leaf = Leaf(self)
        self.stem = Stem(self)

    # @property
    # def initiated(self):
    #     return True

    @property
    def mass(self):
        return self.leaf.mass + self.stem.mass

    #TODO handle predawn_lp elsewhere
    def update(self):
        self.leaf.update()
        self.stem.update()
