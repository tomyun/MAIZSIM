from .trait import Trait

import numpy as np

class Mass(Trait):
    def setup(self):
        # seed weight g/seed
        self._seed = 0.275

    @property
    def seed(self):
        return self._seed

    #TODO handle carbon supply from the seed
    def reduce_seed(self, supply):
        supply = np.fmin(self._seed, supply)
        self._seed -= supply
        return supply

    @property
    def stem(self):
        # dt the addition of C_reserve here only serves to maintain a total for the mass. It could have just as easily been added to total mass.
        # C_reserve is added to stem here to represent soluble TNC, SK
        return sum([nu.stem.mass for nu in self.p.nodal_units]) + self.p.carbon.reserve

    @property
    def initial_leaf(self):
        return self.seed * self.p.ratio.initial_leaf

    # this is the total mass of active leaves that are not entirely dead (e.g., dropped).
    # It would be slightly greather than the green leaf mass because some senesced leaf area is included until they are complely aged (dead), SK
    @property
    def active_leaf(self):
        return sum([nu.leaf.mass for nu in self.p.nodal_units if not nu.leaf.dropped])

    @property
    def dropped_leaf(self):
        return sum([nu.leaf.mass for nu in self.p.nodal_units if nu.leaf.dropped])

    @property
    def total_leaf(self):
        # this should equal to activeLeafMass + droppedLeafMass
        return sum([nu.leaf.mass for nu in self.p.nodal_units])

    @property
    def leaf(self):
        return self.total_leaf

    @property
    def ear(self):
        return self.p.ear.mass

    @property
    def root(self):
        return self.p.root.mass

    @property
    def shoot(self):
        return self.seed + self.stem + self.leaf + self.ear

    @property
    def total(self):
        return self.shoot + self.root

    # this will only be used for total leaf area adjustment.
    # If the individual leaf thing works out this will be deleted.
    @property
    def potential_carbon_demand(self):
        # Just a mocking value for now. Need to find a more mechanistic way to simulate change in SLA YY
        # SK 8/20/10: changed it to 200 cm2/g based on data from Kim et al. (2007) EEB
        SLA = 200

        # units are biomass not carbon
        leaf_mass_demand = self.p.area.potential_leaf_increase / SLA
        # potential_carbon_demand = carbon_demand # for now only carbon demand for leaf is calculated.
        return leaf_mass_demand
