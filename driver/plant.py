from phenology import Phenology
from structure import NodalUnit
from gasexchange import GasExchange

import numpy as np

class PlantTait:
    def __init__(self, plant):
        self.p = plant
        self.setup()

    def setup(self):
        pass


class Ratio(PlantTrait):
    @property
    def carbon_to_mass(self):
        # Weight.C_to_CH2O_ratio
        # 40% C, See Kim et al. (2007) EEB
        return 0.40

    @property
    def shoot_to_root(self):
        return 0.7

    @property
    def root_to_shoot(self):
        return 1 - self.shoot_to_root

    @property
    def leaf_to_stem(self):
        return 0.9

    @property
    def stem_to_leaf(self):
        return 1 - self.leaf_to_stem

    @property
    def initial_leaf(self):
        #TODO how to handle primordia?
        return self.shoot_to_root * self.leaf_to_stem / self.p.primordia


class Mass(PlantTrait):
    def setup(self):
        # seed weight g/seed
        self._seed = 0.275

    @property
    def seed(self):
        return self._seed

    #TODO handle carbon supply from the seed
    @property
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
        return sum([nu.leaf.mass for nu in self.p.nodal_units if not nu.dropped])

    @property
    def dropped_leaf(self):
        return sum([nu.leaf.mass for nu in self.p.nodal_units if nu.dropped])

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


class Area(PlantTrait):
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
    def per_leaf_relative_increase(self):
        return sum([nu.leaf.relative_area_increase for nu in self.p.nodal_units])


class Count(PlantTraint):
    @property
    def total_growing_leaves(self):
        return sum([nu.leaf.growing for nu in self.p.nodal_units])

    @property
    def total_dropped_leaves(self):
        return sum([nu.leaf.dropped for nu in self.p.nodal_units])


class Carbon(PlantTrait):
    def setup(self):
        self.reserve = 0
        self.pool = 0

        self.supply = 0 # daily mobilization of carbon
        self.demand = 0

    def translocate_to_pool(self, amount=None):
        if amount is None:
            amount = self.reserve
        else:
            amount = np.fmin(self.reserve, amount)
        self.reserve -= amount
        self.pool += amount

    def assimilate_to_pool(self, amount=None):
        amount = self.p.photosynthesis.assimilation
        self.pool += amount

    #TODO merge consume_pool / reserve with a solid logic
    def consume_pool(self, amount):
        #FIXME no boundary check
        self.pool -= amount

    def consume_reserve(self, amount):
        self.reserve -= amount

    def reset_pool(self):
        # reset shorterm C_pool to zero at midnight, needs to be more mechanistic
        #FIXME need np.fmax(0, self.pool)?
        self.reserve += self.pool
        self.pool = 0

    def allocate_with_seed(self):
        self.reserve = self.p.mass.seed * self.p.ratio.carbon_to_mass
        # assume it takes 20 days to exhaust seed C reserve
        #self.translocate_to_pool(self.reserve * (1/20) * (1/24) * (initInfo.timeStep / 60))
        self.translocate_to_pool()
        #FIXME the original code did not reset pool here

    # to be used by allocate_carbon()
    @property
    def _temperature_effect(self):
        #FIXME properly handle T_air
        T_air = self.p.atmos.T_air

        # this needs to be f of temperature, source/sink relations, nitrogen, and probably water
        # a valve function is necessary because assimilates from CPool cannot be dumped instantanesly to parts
        # this may be used for implementing feedback inhibition due to high sugar content in the leaves
        # The following is based on Grant (1989) AJ 81:563-571

        # Normalized (0 to 1) temperature response fn parameters, Pasian and Lieth (1990)
        # Lieth and Pasian Scientifica Hortuculturae 46:109-128 1991
        # parameters were fit using rose data -
        b1 = 2.325152587
        b2 = 0.185418876 # I'm using this because it can have broad optimal region unlike beta fn or Arrhenius eqn
        b3 = 0.203535650
        Td = 48.6 #High temperature compensation point

        g1 = 1 + np.exp(b1 - b2 * T_air)
        g2 = 1 - np.exp(-b3 * np.fmax(0, Td - T_air))
        return g2 / g1

    @property
    def _growth_factor(self):
        # translocation limitation and lag, assume it takes 1 hours to complete, 0.2=5 hrs
        # this is where source/sink (supply/demand) valve can come in to play
        # 0.2 is value for hourly interval, Grant (1989)
        hours = 5
        return 1 / (hours * 60 / info.time_step)

    @property
    def _scale(self):
        # see Grant (1989), #of phy elapsed since TI/# of phy between TI and silking
        return self.p.pheno.leaf_appearance_fraction_since_tassel_initiation

    # based on McCree's paradigm, See McCree(1988), Amthor (2000), Goudriaan and van Laar (1994)
    # units very important here, be explicit whether dealing with gC, gCH2O, or gCO2
    @property
    def maintenance_respiration(self):
        Q10 = 2.0 # typical Q10 value for respiration, Loomis and Amthor (1999) Crop Sci 39:1584-1596
        #Q10 = 2.1 # where does this value come from?

        #FIXME proper hanlding of init object
        dt = self.p.init.time_step / (24 * 60)

        # gCH2O g-1DM day-1 at 20C for young plants, Goudriaan and van Laar (1994) Wageningen textbook p 54, 60-61
        #coeff = 0.015
        coeff = 0.018

        # as more leaves senesce maint cost should go down, added 1 to both denom and numer to avoid division by zero.
        #agefn = (self.green_leaf_area + 1) / (self.leaf_area + 1)
        # no maint cost for dead materials but needs to be more mechanistic, SK
        agefn = 1.0

        #FIXME proper handling of atmos object
        q10fn = Q10 ** ((self.p.atmos.air_T - 20.0) / 10) # should be soil temperature
        return q10fn * coeff * self.p.mass.total * dt # gCH2O dt-1, agefn effect removed. 11/17/14. SK.

    #TODO rethink the entire logic!!!
    def make_supply(self):
        # C_demand does not enter into equations until grain fill
        translocation_rate = self._temperature_effect * self._growth_factor
        maintenance_respiration = self.maintenance_respiration

        if self.pool > self.demand:
            # CADD from Grant
            self.supply = self.pool * translocation_rate
            # C_Pool is what is left over from the carbon available for growth
            self.consume_pool(supply)
        elif self.pool == 0: # np.abs(self.pool) < 0.0001:
            if self.reserve > 0:
                # all the reserve is not available
                #FIXME why translocation from the reserve has the same rate?
                self.supply = self.reserve * translocation_rate
                # reduce reserve pool for used carbon
                self.consume_reserve(supply)
            else:
                #TODO what would happen if no reserve available?
                pass
        elif self.reserve > self.demand > 0:
            # conversion and translocation from long term reserve should be less efficient, apply nother glucose conversion factor
            # translocation from the soluble sugar reserve
            #TODO: prevent negative pool beforehand?
            if self.pool < 0:
                # C_pool negative here, add instead of subtract
                self.reset_pool()
            # deplete C_pool first* tmprEffect
            #FIXME why supply size depends on the total demand unlike other cases?
            self.supply = self.demand * translocation_rate
            self.consume_reserve(supply) # reserve C is used to increase grain mass
            #FIXME isn't it duplicate of negative pool handling right above?
            # send remaining C (not enough to meet the demand) in shorterm pool to reserve
            # empty the C_pool
            self.reset_pool()

        #FIXME handling of maintenance respiration starts here, but not sure it should be
        elif self.pool > maintenance_respiration:
            self.supply = maintenance_respiration
            self.consume_pool(maintenance_respiration)
        elif self.reserve > maintenance_respiration:
            self.supply = maintenance_respiration
            self.consume_reserve(maintenance_respiration)
            self.reset_pool()
            # In this way, the c_reserve is used to satisfy the maintainance respiration demand. Then C_pool is
            # used to replenish C_reserve and then C_pool is set to 0. In other words, C_pool is depleted first
            # and then C_reserve is used to supplement whatever demand that's left
        else:
            self.reset_pool()
            self.supply = np.fmin(self.reserve, maintenance_respiration)

    @property
    def partition(self):
        fraction = np.fmin(0.925, 0.50 + 0.50 * self._scale) # eq 3 in Grant
        #conv_factor = 1 / 1.43 # equivalent Yg, Goudriaan and van Laar (1994)
        Yg = 0.75 # synthesis efficiency, ranges between 0.7 to 0.76 for corn, see Loomis and Amthor (1999), Grant (1989), McCree (1988)
        # Yg = 0.74

        c = supply - self.maintenance_respiration
        # this is the same as (PhyllochronsSinceTI - lvsAtTI / (totalLeaves - lvsAtTI))
        return {
            'shoot': np.fmax(0, Yg * fraction * c), # gCH2O partitioned to shoot
            'root': np.fmax(0, Yg * (1 - fraction) * c), # gCH2O partitioned to roots
        }

    @property
    def partition_shoot(self):
        #H
        w.


    def prepare_mobilization(self):




#TODO move into Leaf class?
class Nitrogen(PlantTrait):
    # SK: get N fraction allocated to leaves, this code is just moved from the end of the procedure, this may be taken out to become a separate fn

    def setup(self):
        # assume nitrogen concentration at the beginning is 3.4% of the total weight of the seed
        # need to check with Yang. This doesn't look correct
        self.set_pool(self.p.mass.seed * 0.034)

    def set_pool(self, pool):
        shoot_mass = self.p.mass.shoot
        if shoot_mass * self.p.info.plant_density <= 100: # g m-2
            # when shoot biomass is lower than 100 g/m2, the maximum [N] allowed is 6.3%
            # shoot biomass and Nitrogen are in g
            # need to adjust demand or else there will be mass balance problems
            pool = np.fmin(0.63 * shoot_mass, pool)
        self._pool = pool

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


#TODO rename to CarbonAssimilation or so? could be consistently named as CarbonPartition, CarbonAllocation...
class Photosynthesis(PlantTrait):
    def setup(self):
        self.sunlit_leaf = GasExchange('Sunlit')
        self.shaded_leaf = GasExchange('Shaded')

    @property
    def sunlit_leaf_area_index(self):
        return self.p.lightenv.sunlitLAI()

    @property
    def shaded_leaf_area_index(self):
        return self.p.lightenv.shadedLAI()

    @property
    def leaf_area_index_array(self):
        return np.array([
            self.sunlit_leaf_area_index,
            self.shaded_leaf_area_index,
        ])

    def _weighted(self, array):
        return self.leaf_area_index_array.dot(array)

    @property
    def gross_array(self):
        return np.array([
            self.sunlit.A_gross,
            self.shaded.A_gross,
        ])

    @property
    def net_array(self):
        return np.array([
            self.sunlit.A_net,
            self.shaded.A_net,
        ])

    @property
    def evapotranspiration_array(self):
        return np.array([
            self.sunlit.ET,
            self.shaded.ET,
        ])

    @property
    def temperature_array(self):
        return np.array([
            self.sunlit.T_leaf,
            self.shaded.T_leaf,
        ])

    @property
    def conductance_array(self):
        return np.array([
            self.sunlit.gs,
            self.shaded.gs,
        ])

    @property
    def gross_CO2_umol_per_m2_s(self):
        return self._weighted(self.gross_array)

    # plantsPerMeterSquare units are umol CO2 m-2 ground s-1
    # in the following we convert to g C plant-1 per hour
    # photosynthesis_gross is umol CO2 m-2 leaf s-1

    @property
    def net_CO2_umol_per_m2_s(self):
        # grams CO2 per plant per hour
        return self._weighted(self.net_array)

    @property
    def transpiration_H2O_mol_per_m2_s(self):
        #TODO need to save this?
        # when outputting the previous step transpiration is compared to the current step's water uptake
        #self.transpiration_old = self.transpiration
        #FIXME need to check if LAIs are negative?
        #transpiration = sunlit.ET * max(0, sunlit_LAI) + shaded.ET * max(0, shaded_LAI)
        return self._weighted(self.evapotranspiration_array)

    #TODO consolidate unit conversions somewhere else

    @property
    def _mol_per_umol(self):
        return 1 / 1e6

    @property
    def _plant_per_m2(self):
        return 1 / initInfo.plant_density

    @property
    def _min_step_per_sec(self):
        return 60 * initInfo.time_step

    # final values

    @property
    def assimilation(self):
        # grams CO2 per plant per hour
        return np.prod([
            self.gross_CO2_umol_per_m2_s,
            self._mol_per_umol,
            self._plant_per_m2,
            self._min_step_per_sec,
            Weight.CO2,
        ])

    @property
    def gross(self):
        # grams carbo per plant per hour
        return np.prod([
            self.gross_CO2_umol_per_m2_s,
            self._mol_per_umol,
            self._plant_per_m2,
            self._min_step_per_sec,
            Weight.CH2O,
        ])

    @property
    def net(self):
        # grams carbo per plant per hour
        return np.prod([
            self.net_CO2_umol_per_m2_s,
            self._mol_per_umol,
            self._plant_per_m2,
            self._min_step_per_sec,
            Weight.CH2O,
        ])

    @property
    def transpiration(self):
        # Units of Transpiration from sunlit->ET are mol m-2 (leaf area) s-1
        # Calculation of transpiration from ET involves the conversion to gr per plant per hour
        #FIXME _min_step_per_sec used instead of fixed 3600 = 60 * 60
        return np.prod([
            self.transpiration_H2O_mol_per_m2_s,
            self._plant_per_m2,
            self._min_step_per_sec,
            Weight.H2O,
        ])

    @property
    def temperature(self):
        return self._weighted(self.temperature_array)

    #TODO is it needed?
    @property
    def vapor_pressure_deficit(self):
        return self.sunlit.VPD

    @property
    def conductance(self):
        #TODO is this condition necessary?
        #if sunlit_LAI >= 0 and shaded_LAI >= 0 and LAI >= 0:
        try:
            # average stomatal conductance Yang
            c = self._weighted(self.conductance_array) / self.p.area.leaf_area_index
            return np.fmax(0, c)
        except ZeroDivisionError:
            return 0


#TODO split into multiple mixins
class Plant:
    def __init__(self, info):
        #TODO implement InitInfo alternatives
        self.info = info

        #timestep = info...
        #TODO pass PRIMORDIA as initial_leaves
        self.primordia = 5

        self.pheno = Phenology(timestep)
        self.mass = Mass(self)
        self.area = Area(self)
        self.count = Count(self)
        self.carbon = Carbon(self)
        self.nitrogen = Nitrogen(self)
        self.photosynthesis = Photosynthesis(self)

        #TODO make another trait object for the structure?
        self.setup_structure()

    def setup_structure(self):
        self.root = None
        self.ear = None
        self.nodal_units = []

    def initiate_primordia(self):
        self.nodal_units = [
            NodalUnit(
                self,
                rank=i+1,
                mass=self.mass.initial_leaf
            ) for i in range(self.primordia)
        ]

    def initiate_root(self):
        # here we calculate the mass of roots that were initialized in the soil model (read in with the element data)
        # This is so there is no discontinuity in root mass (relative to the carbon supplied by the plant later)
        # these are roots taht grew from the seed

        self.root = Root()
        #TODO use weather.TotalRootWeight from 2DSOIL
        self.root.import_carbohydrate(soil.total_root_weight)

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

    def update(self, atmos, predawn_lwp):
        #TODO pass atmos as is?
        self.pheno.update(atmos.T_air)

        if pheno.germinating():
            #TODO temperature setting?
        elif pheno.emerging():
            self.initiate_root()

            #TODO temperature setting?

            #TODO carbon reserve/pool init

            #HACK update_leaves() should precede initiate_leaves()
            self.update_leaves()
            self.initiate_leaves()

            # calculate relative area increases for leaves now that they are updated
            # commented for now, have to test this for seedling
            #calcMaintRespiration(weather)

            #TODO logics clean up
            if not first_leaf.appeared:
                self.maintenance_respiration(atmos) #FIXME no side-effect
                self.allocate_carbon()
                self.update_seed_mass()
            else:
                self.calc_gas_exchange()
                self.carbon.assimilate_to_pool()
                self.maintenance_respiration(atmos) #FIXME no side-effect
                self.allocate_carbon()
                self.update_seed_mass() #FIXME with different ratio
        elif not pheno.dead():
            self.calc_gas_exchange()
            self.carbon.assimilate_to_pool()
            #self.maintenance_respiration(atmos) #FIXME no side-effect
            self.allocate_carbon()
            #TODO need DateTime like object
            if atmos.time == midnight:
                self.carbon.reset_pool()
            else:
                self.carbon.assimilate_to_pool()
            #self.update_mass()


    ###########
    # Process #
    ###########

    #TODO unify naming convention for processes (i.e. verb?)

    def calc_gas_exchange(self, atmos):
        #tau = 0.50 # atmospheric transmittance, to be implemented as a variable => done

        atm_pressure = 100.0 # kPa, to be predicted using altitude
        #FIXME fix atmos.P_air?
        atmos.P_air = atm_pressure

        LAF = 1.37 # leaf angle factor for corn leaves, Campbell and Norman (1998)
        leaf_width = 5.0 # to be calculated when implemented for individal leaves

        #TODO handle plant_density from initInfo object
        LAI = self.area.leaf_area_index

        #TODO how do we get LeafWP and ET_supply?
        LWP = weather.LeafWP
        ET_supply = weather.ET_supply * initInfo.plant_density / 3600 / 18.01 / LAI

        #TODO integrate lightenv with Atmosphere class?
        #TODO lightenv.dll needs to be translated to C++. It slows down the execution, 3/16/05, SK
        self.lightenv.radTrans2(weather.jday, weather.time, initInfo.latitude, initInfo.longitude, weather.solRad, weather.PFD, LAI, LAF)
        # temp7 = lightenv.getNIRtot()

        # Calculating transpiration and photosynthesis without stomatal control Y
        # call SetVal_NC()

        # Calculating transpiration and photosynthesis with stomatal controlled by leaf water potential LeafWP Y
        self.sunlit.set_val_psil(
             lightenv.sunlitPFD(),
             atmos.T_air, atmos.CO2, atmos.RH, atmos.wind, atmos.P_air,
             self.nitrogen.leaf_content, leaf_width, LWP, ET_supply
        )

        self.shaded.set_val_psil(
             lightenv.shadedPFD(),
             atmos.T_air, atmos.CO2, atmos.RH, atmos.wind, atmos.P_air,
             self.nitrogen.leaf_content, leaf_width, LWP, ET_supply
        )

    def allocate_carbon(self, atmos):
        supply = self.carbon.supply()
