from ..phenology import Phenology
from ..morphology import Ear, Root, NodalUnit
from .mass import Mass
from .area import Area
from .count import Count
from .ratio import Ratio
from .carbon import Carbon
from .nitrogen import Nitrogen
from .water import Water
from .photosynthesis import Photosynthesis
from ..timer import Timer

import numpy as np

#TODO split into multiple mixins
class Plant:
    def __init__(self, initials, variety):
        self.initials = initials
        self.variety = variety

        #TODO pass PRIMORDIA as initial_leaves
        self.primordia = 5

        self.weather = None
        self.soil = None

        # phenology
        self.pheno = Phenology(self)

        # morphology
        #TODO make another trait object for the structure?
        self.setup_structure()

        # physiological traits
        self.mass = Mass(self)
        self.area = Area(self)
        self.count = Count(self)
        self.ratio = Ratio(self)
        self.carbon = Carbon(self)
        self.nitrogen = Nitrogen(self)
        self.water = Water(self)
        self.photosynthesis = Photosynthesis(self)

        #HACK initiate primordia here
        self.initiate_primordia()

    def setup_structure(self):
        #TODO initiate_ear() like leaves and root
        self.ear = Ear(self)
        self.root = Root(self)
        self.nodal_units = []

    def initiate_primordia(self):
        self.nodal_units = [NodalUnit(self, rank=i+1) for i in range(self.primordia)]
        #FIXME should not allocate leaf mass here, otherwise accounting total seed mass gets complicated
        #[nu.leaf.import_carbohydrate(self.mass.initial_leaf) for nu in self.nodal_units]

    def initiate_root(self):
        # here we calculate the mass of roots that were initialized in the soil model (read in with the element data)
        # This is so there is no discontinuity in root mass (relative to the carbon supplied by the plant later)
        # these are roots taht grew from the seed

        if not self.root.initiated:
            self.root.initiated = True
            self.root.import_carbohydrate(self.soil.total_root_weight)

    def initiate_leaves(self):
        for i in range(self.pheno.leaves_initiated):
            try:
                self.nodal_units[i]
            except IndexError:
                nu = NodalUnit(self, rank=i+1)
                self.nodal_units.append(nu)

    def update_leaves(self):
        #SK 8/22/10: set leaf N content before each nodal unit is updated

        # Pass the predawn leaf water potential into a nodel
        # to enable the model to simulate leaf expansion with the
        # effect of predawn leaf water potential YY

        [nu.update() for nu in self.nodal_units]

    ##########
    # Update #
    ##########

    def update(self, weather, soil):
        self.weather = weather
        self.soil = soil

        self.pheno.update()

        if self.pheno.germinating:
            #TODO temperature setting?
            pass
        elif self.pheno.emerging:
            self.initiate_root()

            #TODO temperature setting?

            #HACK carbon reserve/pool init from seed done in Carbon object

            #HACK update_leaves() should precede initiate_leaves()
            self.update_leaves()
            self.initiate_leaves()

            # calculate relative area increases for leaves now that they are updated
            # commented for now, have to test this for seedling
            #calcMaintRespiration(weather)

            #TODO logics clean up
            first_leaf = self.nodal_units[0].leaf
            if not first_leaf.appeared:
                #self.maintenance_respiration() #FIXME no side-effect
                self.allocate_carbon()
                #self.update_seed_mass()
            else:
                self.calc_gas_exchange()
                self.carbon.assimilate_to_pool()
                #self.maintenance_respiration() #FIXME no side-effect
                self.allocate_carbon()
                #self.update_seed_mass() #FIXME with different ratio
        elif not self.pheno.dead:
            #HACK update_leaves() should precede initiate_leaves()
            self.update_leaves()
            self.initiate_leaves()

            self.calc_gas_exchange()
            self.carbon.assimilate_to_pool()
            #self.maintenance_respiration() #FIXME no side-effect
            self.allocate_carbon()
            if weather.time.hour == 0: # midnight
                self.carbon.reset_pool()
            else:
                self.carbon.assimilate_to_pool()
            #self.update_mass()


    ###########
    # Process #
    ###########

    #TODO unify naming convention for processes (i.e. verb?)

    def calc_gas_exchange(self):
        self.photosynthesis.update()

    def allocate_carbon(self):
        self.carbon.make_supply()

        #FIXME just for validating total imports of carbohydrate
        total_demand = 0

        for i in range(self.pheno.leaves_initiated):
            nu = self.nodal_units[i]

            # here we can allocate leaf part among leaves
            # need to find demand first
            leaf = nu.leaf

            #HACK setting SLA this way is not permitted in the new implementation; this was unused code anyways
            # # Update SLA based on current shoot growth rate. Do this for every leaf
            # # no sla adjustment until after emergence
            # if self.carbon.shoot > 0 and self.pheno.emerging:
            #     # see Grant, 1989, eq 9. 10000 converts from m2 to cm2
            #     SLA_est = 1 / (25.0 + 150.0 * self.carbon.shoot) * 10000
            #     if leaf.growing:
            #         leaf.specific_leaf_area = np.fmin(400.0, SLA_est)

            # Partition carbon to leaves relative to growth rate
            # now need to find increment of Carbo to add to each leaf
            # update leaf's mass
            # first find the number of growing leaves. This is saved as a variable in the [0] nodal unit

            if not leaf.dead:
                # Adjusting C allocation based on leaf size if not aging.
                # doing it based on current growth rate is more mechanistic but seems to have issue now. To revisit. SK
                demand = leaf.potential_area / self.area.potential_leaf * self.carbon.leaf
                # carbon allocated according to growth rate
                #demand = leaf.relative_area_increase / self.area.relative_leaf_increase * self.carbon.leaf
                total_demand += demand

                leaf.import_carbohydrate(demand)
        #assert total_demand == self.carbon.leaf

        #FIXME what is the difference between import_carbohydrate()?
        #self.root.actual_carbon_increment = self.carbon.root
        # before emergence root weight has been initialized. Just dump this carbon for now.
        if self.pheno.emerged:
            self.root.import_carbohydrate(self.carbon.root)

        #FIXME in the original code, total stem mass was stored in the 0th nodal unit
        #TODO partition into individual stem like we did for the leaf
        self.nodal_units[0].stem.import_carbohydrate(self.carbon.stem)

        self.ear.import_carbohydrate(self.carbon.ear)

        # checking the balance if sums up to shootPart
        #part_sum = self.carbon.stem + self.carbon.ear + self.carbon.leaf
        #if self.weather.time.hour == 12: # noon?
            #print("adding CH2O: {} to grain".format(self.carbon.grain))
            #print("Sum of part coeff is {} and shoot CH2O is {}".format(part_sum, self.carbon.shoot))
        #print("Carbon pool = {}".format(self.carbon.pool))
