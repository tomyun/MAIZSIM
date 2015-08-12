from phenology import Phenology
from structure import NodalUnit
from gasexchange import GasExchange

import numpy as np

class PlantTait:
    def __init__(self, plant):
        self.p = plant


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
    @property
    def seed(self):
        # seed weight g/seed
        return 0.275

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
    def active_leaf(self):
        return self.green_leaf / self.leaf

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
    @property
    def reserve(self):
        if self.p.pheno.emerging:
            return self.p.mass.seed * self.p.ratio.carbon_to_mass
        pass

    @property
    def pool(self):
        if self.p.pheno.emerging:
            # assume it takes 20 days to exhaust seed C reserve
            # C_pool += C_reserve * (1/20) * (1/24) * (initInfo.timeStep / 60)
            return self.reserve


class Nitrogen(PlantTrait):
    # SK: get N fraction allocated to leaves, this code is just moved from the end of the procedure, this may be taken out to become a separate fn

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
        return np.max(0, fraction)

    @property
    def leaf(self):
        # calculate total nitrogen amount in the leaves YY units are grams N in all the leaves
        return self.leaf_fraction * self.p.nitrogen

    @property
    def leaf_content(self):
        #SK 8/22/10: set avg greenleaf N content before update in g/m2
        return self.leaf / (self.p.area.green_leaf / 10000)


#TODO split into multiple mixins
class Plant:
    def __init__(self, info):
        #timestep = info...
        #TODO pass PRIMORDIA as initial_leaves
        self.primordia = 5

        self.pheno = Phenology(timestep)
        self.mass = Mass(self)
        self.area = Area(self)
        self.count = Count(self)
        self.carbon = Carbon(self)

        self.setup_structure()

    def setup_structure(self):
        self.root = None
        self.ear = None
        self.nodal_units = []

    def setup_nitrogen(self):
        # assume nitrogen concentration at the beginning is 3.4% of the total weight of the seed
        # need to check with Yang. This doesn't look correct
        self.nitrogen = self.mass.seed * 0.034

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

    ###########
    # Process #
    ###########

    def calc_gas_exchange(self, atmos):
        #tau = 0.50 # atmospheric transmittance, to be implemented as a variable => done

        atm_pressure = 100.0 # kPa, to be predicted using altitude
        #FIXME fix atmos.P_air?
        atmos.P_air = atm_pressure

        LAF = 1.37 # leaf angle factor for corn leaves, Campbell and Norman (1998)
        leaf_width = 5.0 # to be calculated when implemented for individal leaves

        #TODO handle plant_density from initInfo object
        LAI = self.area.green_leaf * initInfo.plant_density / 100**2

        #TODO how do we get LeafWP and ET_supply?
        LWP = weather.LeafWP
        ET_supply = weather.ET_supply * initInfo.plant_density / 3600 / 18.01 / LAI

        sunlit = GasExchange('Sunlit', self.nitrogen.leaf_content)
        shaded = GasExchange('Shaded', self.nitrogen.leaf_content)

        #TODO lightenv.dll needs to be translated to C++. It slows down the execution, 3/16/05, SK
        lightenv.radTrans2(weather.jday, weather.time, initInfo.latitude, initInfo.longitude, weather.solRad, weather.PFD, LAI, LAF)
        # temp7 = lightenv.getNIRtot()

        # Calculating transpiration and photosynthesis without stomatal control Y
        # call SetVal_NC()

        # Calculating transpiration and photosynthesis with stomatal controlled by leaf water potential LeafWP Y
        sunlit.set_val_psil(
             lightenv.sunlitPFD(),
             atmos.T_air, atmos.CO2, atmos.RH, atmos.wind, atmos.P_air,
             leaf_width, LWP, ET_supply
        )

        shaded.set_val_psil(
             lightenv.shadedPFD(),
             atmos.T_air, atmos.CO2, atmos.RH, atmos.wind, atmos.P_air,
             leaf_width, LWP, ET_supply
        )

        sunlit_LAI = lightenv.sunlitLAI()
        shaded_LAI = lightenv.shadedLAI()

        #TODO need to save this?
        # plantsPerMeterSquare units are umol CO2 m-2 ground s-1
        self.photosynthesis_gross = sunlit.A_gross * sunlit_LAI + shaded.A_gross * shaded_LAI
        self.photosynthesis_net = sunlit.A_net * sunlit_LAI + shaded.A_net * shaded_LAI

        # photosynthesis_gross is umol CO2 m-2 leaf s-1
        # in the following we convert to g C plant-1 per hour
        unit = (60.0 * initInfo.time_step) / initInfo.plant_density / 1.0e6
        self.assimilate = photosynthesis_gross * Weight.CO2 * unit # grams CO2 per plant per hour
        self.photosynthesis_gross *= Weight.CH2O * unit # grams carbo per plant per hour
        self.photosynthesis_net *= Weight.CH2O * unit # grams carbo per plant per hour

        #TODO need to save this?
        # when outputting the previous step transpiration is compared to the current step's water uptake
        self.transpiration_old = self.transpiration
        #FIXME need to check if LAIs are negative?
        transpiration = sunlit.ET * max(0, sunlit_LAI) + shaded.ET * max(0, shaded_LAI)
        # plantsPerMeterSquare units are grams per plant per hour
        transpiration /= initInfo.plant_density * 3600 * 18.01

        #TODO need to save this?
        self.VPD = sunlit.VPD

        #TODO is this condition necessary?
        #if sunlit_LAI >= 0 and shaded_LAI >= 0 and LAI >= 0:
        try:
            # average stomatal conductance Yang
            self.conductance = np.max(0, (sunlit.gs * sunlit_LAI + shaded.gs * shaded_LAI) / LAI)
        except ZeroDivisionError:
            self.conductance = 0

    def carbon_allocation(self, atmos):
        pass

    # based on McCree's paradigm, See McCree(1988), Amthor (2000), Goudriaan and van Laar (1994)
    # units very important here, be explicit whether dealing with gC, gCH2O, or gCO2
    def maintenance_respiration(self, atmos):
        Q10 = 2.0 # typical Q10 value for respiration, Loomis and Amthor (1999) Crop Sci 39:1584-1596
        #Q10 = 2.1 # where does this value come from?

        dt = initInfo.time_step / (24 * 60)

        # gCH2O g-1DM day-1 at 20C for young plants, Goudriaan and van Laar (1994) Wageningen textbook p 54, 60-61
        #coeff = 0.015
        coeff = 0.018

        # as more leaves senesce maint cost should go down, added 1 to both denom and numer to avoid division by zero.
        #agefn = (self.green_leaf_area + 1) / (self.leaf_area + 1)
        # no maint cost for dead materials but needs to be more mechanistic, SK
        agefn = 1.0

        q10fn = Q10 ** ((atmos.air_T - 20.0) / 10) # should be soil temperature
        return q10fn * coeff * self.mass.total * dt # gCH2O dt-1, agefn effect removed. 11/17/14. SK.
