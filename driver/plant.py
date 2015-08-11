from phenology import Phenology
from structure import NodalUnit

import numpy as np

class Mass(object):
    def __init__(self, plant):
        self.plant = plant

    @property
    def seed(self):
        # seed weight g/seed
        return 0.275

    @property
    def stem(self):
        pass

    @property
    def leaf(self):
        pass

    @property
    def initial_leaf(self):
        return self.seed * self.plant.initial_leaf_ratio

    @property
    def active_leaf(self):
        pass

    @property
    def dropped_leaf(self):
        pass

    @property
    def ear(self):
        pass

    @property
    def root(self):
        pass

    @property
    def shoot(self):
        return self.seed + self.stem + self.leaf + self.ear

    @property
    def total(self):
        return self.shoot + self.root


class Plant(object):
    def __init__(self, info):
        #timestep = info...
        #TODO pass PRIMORDIA as initial_leaves
        self.primordia = 5
        self.pheno = Phenology(timestep)
        self.setup_structure()

    def setup_structure(self):
        self.roots = Roots()
        self.ear = Ear()
        self.nodal_units = [
            NodalUnit(
                self.pheno,
                rank=i,
                mass=self.mass.initial_leaf
            ) for i in range(self.primordia)
        ]

    # Mass

    @property
    def shoot_to_root_ratio(self):
        return 0.7

    @property
    def root_to_shoot_ratio(self):
        return 1 - self.shoot_ratio

    @property
    def leaf_to_stem_ratio(self):
        return 0.9

    @property
    def stem_to_leaf_ratio(self):
        return 1 - self.leaf_to_stem_ratio

    @property
    def initial_leaf_ratio(self):
        return self.shoot_to_root_ratio * self.leaf_to_stem_ratio / self.primordia

    def update(self, atmos, predawn_lwp):
