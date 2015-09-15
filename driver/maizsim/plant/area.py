from .trait import Trait

class Area(Trait):
    @property
    def leaf(self):
        return sum([nu.leaf.area for nu in self.p.nodal_units])

    @property
    def green_leaf(self):
        return sum([nu.leaf.green_area for nu in self.p.nodal_units])

    #TODO remove if unnecessary
    @property
    def active_leaf_ratio(self):
        return self.green_leaf / self.leaf

    @property
    def leaf_area_index(self):
        #TODO handle info.plant_density
        return self.green_leaf * info.plant_density / 100**2

    # actualgreenArea is the green area of leaf growing under carbon limitation
	#SK 8/22/10: There appears to be no distinction between these two variables in the code.
    @property
    def actual_green_leaf(self):
        return self.green_leaf

    @property
    def senescent_leaf(self):
        return sum([nu.leaf.senescent_area for nu in self.p.nodal_units])

    @property
    def potential_leaf(self):
        return sum([nu.leaf.potential_area for nu in self.p.nodal_units])

    @property
    def potential_leaf_increase(self):
        return sum([nu.leaf.potential_area_increase for nu in self.p.nodal_units])

    # calculate relative area increases for leaves now that they are updated
    @property
    def relative_leaf_increase(self):
        return sum([nu.leaf.relative_area_increase for nu in self.p.nodal_units])
