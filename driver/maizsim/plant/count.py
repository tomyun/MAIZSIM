from .trait import Trait

class Count(Traint):
    @property
    def total_growing_leaves(self):
        return sum([nu.leaf.growing for nu in self.p.nodal_units])

    @property
    def total_dropped_leaves(self):
        return sum([nu.leaf.dropped for nu in self.p.nodal_units])
