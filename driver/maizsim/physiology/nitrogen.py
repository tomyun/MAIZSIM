from .trait import Trait

import numpy as np

#TODO move into Leaf class?
class Nitrogen(Trait):
    # SK: get N fraction allocated to leaves, this code is just moved from the end of the procedure, this may be taken out to become a separate fn

    def setup(self):
        # assume nitrogen concentration at the beginning is 3.4% of the total weight of the seed
        # need to check with Yang. This doesn't look correct
        self.set_pool(self.initial_pool)

        #TODO set up interface
        self.ratio = 0
        self.hourly_soil_uptake = 0
        self.hourly_demand = 0
        self.cumulative_demand = 0
        self.cumulative_soil_uptake = 0

    def set_pool(self, pool):
        shoot_mass = self.p.mass.shoot
        if shoot_mass * self.p.initials.plant_density <= 100: # g m-2
            # when shoot biomass is lower than 100 g/m2, the maximum [N] allowed is 6.3%
            # shoot biomass and Nitrogen are in g
            # need to adjust demand or else there will be mass balance problems
            #FIXME self.initial_pool or just pool?
            pool = np.fmin(0.063 * shoot_mass, self.initial_pool)
        self._pool = pool

    @property
    def initial_pool(self):
        # assume nitrogen concentration at the beginning is 3.4% of the total weight of the seed
        return 0.034 * self.p.mass.initial_seed # 0.275

    @property
    def pool(self):
        return self._pool

    #TODO for 2DSOIL interface
    def uptake_from_soil(self, amount):
        self.set_pool(self.pool + amount)

    #TODO currently not implemented in the original code
    def remobilize(self):
        pass
        #droppedLfArea = (1-greenLeafArea/potentialLeafArea)*potentialLeafArea; //calculated dropped leaf area YY
        #SK 8/20/10: Changed it to get all non-green leaf area
        #currentDroppedLfArea = droppedLfArea - previousDroppedlfArea; //leaf dropped at this time step
        #this->set_N((this->get_N()-(leaf_N/leafArea)*currentDroppedLfArea)); //calculated the total amount of nitrogen after the dropped leaves take some nitrogen out
        #no nitrogen remobilization from old leaf to young leaf is considered for now YY

    #TODO rename to `leaf_to_plant_ratio`? or just keep it?
    @property
    def leaf_fraction(self):
        # Calculate faction of nitrogen in leaves (leaf NFraction) as a function of thermal time from emergence
        # Equation from Lindquist et al. 2007 YY
        #SK 08/20/10: TotalNitrogen doesn't seem to be updated at all anywhere else since initialized from the seed content
        #SK: I see this is set in crop.cpp ln 253 from NUptake from 2dsoil
        # but this appears to be the amount gained from the soil for the time step; so how does it represent totalNitrogen of a plant?

        # record thermal time from emergency YY
        tt = self.p.pheno.gdd_after_emergence

        # Calculate faction of nitrogen in leaves (leaf NFraction) as a function of thermal time from emergence
        # Equation from Lindquist et al. 2007 YY
        fraction = 0.79688 - 0.00023747 * tt - 0.000000086145 * tt**2

        # fraction of leaf n in total shoot n can't be smaller than zero. YY
        return np.fmax(0, fraction)

    #TODO rename to `leaves`?
    @property
    def leaf(self):
        # calculate total nitrogen amount in the leaves YY units are grams N in all the leaves
        return self.leaf_fraction * self.pool

    #TODO rename to `unit_leaf`?
    # Calculate leaf nitrogen content of per unit area
    @property
    def leaf_content(self):
        # defining leaf nitrogen content this way, we did not consider the difference in leaf nitrogen content
        # of sunlit and shaded leaf yet YY
        #SK 8/22/10: set avg greenleaf N content before update in g/m2
        return self.leaf / (self.p.area.green_leaf / (100**2))
