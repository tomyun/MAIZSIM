import numpy as np

import crop
import soil

# shortcuts for 2DSOIL common blocks
S = soil.shootr
W = soil.weath
G = soil.grid_public
N = soil.nodal_public
#E = soil.elem_public
#B = soil.bound_public
T = soil.time_public
M = soil.module_public
F = soil.datafilenames

from maizsim.io import config
from maizsim.controller import Controller
from maizsim.atmosphere import Weather
from maizsim.rhizosphere import Soil

class Driver:
    def __init__(self, runfile='test.dat'):
        #FIXME [:] no longer works with Python 3?
        soil.datafilenames.runfile[:len(runfile)] = runfile
        self.setup(runfile)

    def setup(self, runfile):
        self.load_config(runfile)

    #TODO implement more consolidated management
    def load_config(self, runfile):
        r = config.Run(runfile)
        self.initials = config.Initials(r.initials)
        self.variety = config.Variety(r.variety)
        self.time = config.Time(r.time)

    def _init(self):
        soil.initialize()
        soil.get_grid_and_boundary()

    def _preprocess(self):
        soil.synchronizer()

        if sum([bool(f) for f in [T.hourlyweather, T.dailyweather]]) != 1:
            raise Exception("error in weather file type")

        if T.hourlyweather:
            soil.setsurfaceh()
        if T.dailyweather:
            soil.setsurfaced()

        soil.settdb()
        soil.autoirrigate()
        soil.mngm()

    def _process(self):
        F.starter = 1.

        if T.linput == 1:
            self._setup_controller(F, S, W, T)
            self._setup_shootr(S)
            self._reset_uptakes(S)
            self._setup_time(T, M)
            self._setup_vars()

        if M.nshoot > 0:
            self._update_uptakes(S, T)

        if abs(T.time - T.tnext[self.mod_num - 1]) < abs(0.001 * T.step):
            # If the sowing date has come and there is not plant, let the program know so other calculations are not done.
            if M.nshoot == 0 and abs(T.time - self.initials.sowing_day) < 0.001:
                M.nshoot = 1

            self._update_nitrogen_uptake_error(S)

            #HACK replace TWeather with Weather, Soil objects
            #w = self._create_weather(T, W)
            #self._update_weather_waterpotential(w, S)
            #self._update_weather_root(w, S, N)
            #self._update_weather_evapotranspiration(w, S)
            #self._update_weather_soil(w, G, N)
            #FIXME only remaining is ET_supply update
            self._update_plant_evapotranspiration(S)

            #self._run_controller(w, self.predawn_lwp)
            self._run_controller()

            if not self.plant.pheno.germinated:
                self._handle_not_germinated(w)

            if self.plant.pheno.emerged:
                self._handle_emerged(w, S)

            self._handle_dead_or_not(M, T)

    def _setup_controller(self, F, S, W, T):
        # Parse the file names from the FORTRAN strings passed from 2dsoil.
        def get_filenames(F):
            vf = F.varietyfile.tostring().strip()
            gf = F.plantgraphics.tostring().strip()
            lf = F.leafgraphics.tostring().strip()
            return vf, gf, lf

        # Set up parameters from initials file whose reading is moved to soil model (Init.for).
        # def create_init_info(S, W, T):
        #     ii = crop.TInitInfo()
        #     ii.plantDensity = S.poparea.item()
        #     ii.latitude = W.latude.item()
        #     ii.longitude = W.longitude.item()
        #     ii.altitude = W.altitude.item()
        #     ii.year = T.year.item()
        #     ii.sowingDay = T.sowingday.item()
        #     ii.beginDay = T.beginday.item()
        #     ii.endDay = T.endday.item()
        #     ii.timeStep = T.timestep.item()
        #     return ii

        vf, gf, lf = get_filenames(F)
        #ii = create_init_info(S, W, T)
        #self.controller = crop.CController(vf, gf, lf, ii)
        #TODO integrate load_config with Controller
        self.controller = Controller(self.initials, self.variety, self.time, gf, lf)

    def _setup_shootr(self, S):
        S.lcai = 0
        S.lareat = 0
        S.height = 0
        # change for debugging purposes
        S.convr = 1 # was 0.1 or 0.38 should be 1.0 as dry matter is used in all cases
        S.awups = 0 # initialize AWUPS, AWUPS_old and LeafWP in 2DSOIL
        S.leafwp = -0.5
        S.pcrs = 0
        S.et_demand = 0
        S.hourlycarboused = 0 # it is also zero'd upon initialization in 2dsoil

    @property
    def period(self):
        # period should be in days, input in minutes
        return T.timestep / 60. / 24.

    @property
    def pop_slab(self):
        #plant_density = S.poprow * 100. / S.rowsp
        #return S.poprow / 100. * S.rowsp * S.eomult
        return S.poprow / 100. * S.eomult

    @property
    def plant(self):
        return self.controller.plant

    # @property
    # def develop(self):
    #     return self.plant.get_develop()

    # @property
    # def initinfo(self):
    #     return self.controller.getInitInfo()

    #########
    # Misc. #
    #########

    #FIXME: is this a proper name?
    def _setup_time(self, T, M):
        T.itime = 1

        T.runflag = 1

        # module number
        M.nummod += 1
        self.mod_num = M.nummod.item()

        T.tnext[self.mod_num - 1] = self.initials.sowing_day

    #FIXME: can we eliminate them?
    def _setup_vars(self):
        # SLNmin: base Specific leaf nitrogen content; for now assume it's 0.5 YY
        self.sln_min = 0.5
        self.predawn_lwp = -0.05
        self.old_shoot_weight_per_m2 = 0.
        self.nitrogen_uptake_old = 0.
        self.cumulative_nitrogen_demand = 0. # grams plant-1

    ###########
    # Uptakes #
    ###########

    def _reset_uptakes(self, S):
        # nitrogen uptake value from 2dsoil accumulated between time steps mg/plant
        #SK 8/20/10: This is curious but OK
        self.nitrogen_uptake = self.plant.nitrogen.pool * self.pop_slab

        # hourly water uptake from 2dsoil
        self.water_uptake = 0.

        S.ndemanderror = self.nitrogen_demand_error = 0
        S.cumulativendemanderror = self.cumulative_nitrogen_demand_error = 0

    def _update_uptakes(self, S, T):
        self._update_nitrogen_uptake(S)
        self._update_water_uptake(S, T)

    def _update_nitrogen_uptake(self, S):
        # Note that SIncrSink has been multiplied by time step in the solute uptake routing the 1000 scales from ug to mg.
        self.nitrogen_uptake += S.sincrsink / 1e6

        if self.nitrogen_uptake > 0:
            #FIXME currently no handling of nitrogen loss
            nuptake = self.nitrogen_uptake / self.pop_slab
            leafloss = self.plant.area.dropped_leaf
            nloss = leafloss * self.sln_min

            #SK 8/20/10: Here seems to be the only place where totalN of the plant is set.
            # NitrogenUptake is initiated from get_N at the begining of the timestep so OK.
            self.plant.nitrogen.set_pool(nuptake)
            # Units are converted from g slab-1 to g plant -1 YY
            # need to look at loss of N in the dropped leaf (plant N goes negative?)
            #FIXME no nitrogen handling in Plant
            #self.plant.set_N(nuptake - nloss)

    def _update_water_uptake(self, S, T):
        self.water_uptake += S.awups * T.step

    def _update_nitrogen_uptake_error(self, S):
        # Calculate error for demand and actual uptake, if negative, demand is greater then uptake.
        #FIXME no nitrogen handling in Plant
        #err = self.nitrogen_uptake / self.pop_slab - self.plant.get_CumulativeNitrogenDemand()
        err = self.nitrogen_uptake / self.pop_slab - self.cumulative_nitrogen_demand_error
        self.nitrogen_demand_error = err
        self.cumulative_nitrogen_demand_error += err
        S.ndemanderror = self.nitrogen_demand_error
        S.cumulativendemanderror = self.cumulative_nitrogen_demand_error

    ############
    # TWeather #
    ############

    def _create_weather(self, T, W):
        i = T.itime - 1
        w = crop.TWeather()

        w.HourlyOutput = T.hourlyoutput.item()
        w.DailyOutput = T.dailyoutput.item()
        w.jday = W.jday.item()
        w.time = T.time - W.jday
        w.CO2 = W.co2.item()
        w.airT = W.tair[i].item()
        w.PFD = W.par[i] * 4.6 # conversion from PAR in Wm-2 to umol s-1 m-2
        w.solRad = W.wattsm[i].item() # conversion from Wm-2 to J m-2 in one hour Total Radiation incident at soil surface

        Es = 0.611 * np.exp(17.502 * w.airT / (240.97 + w.airT)) # saturated vapor pressure at airT
        w.RH = (1 - W.vpd[i] / Es) * 100.

        w.rain = W.rint[i].item()
        w.wind = W.wind * (1000. / 3600.) # conversion from km hr-1 to m s-1
        w.dayLength = W.daylng.item()

        return w

    def _update_weather_waterpotential(self, w, S):
        # since LeafWP in 2dsoil is in bar but in maizesim is in MPa, so, have to divide it by 10 to convert it into MPa before passing the value to Maizesim 1 bar=10kPa
        w.LeafWP = S.psil_ / 10. # and leaf water potential information into MAIZESIM Yang 8/15/06

        # If time is 5 am, then pass the leaf water potential (the predawn leaf water potential) from SHOOTR to the wthr object. YY
        if abs(w.time - 0.2083) < 0.0001:
            # Here LeafWP is in bar.
            # Since the LWPeffect in leaf.cpp uses leaf water potential in bar,
            # so here PredawnLWP is in bar, instead of being scaled to MPa. YY
            self.predawn_lwp = S.psil_.item() # SHOOTR->LeafWP

    def _update_weather_root(self, w, S, N):
        w.pcrl = S.pcrl / self.pop_slab / 24.
        w.pcrq = S.pcrq / self.pop_slab / 24.

        # Pass actual carbohydrate amount used in 2dsoil back to the plant.
        #ToDo - make pcrs a new variable (ActualRootCarboUsed) and make it a member of plant.
        #dt here I changed this temporarily for debugging
        # don't need to divide by 24 since the value has been integrated over an hour
        w.pcrs = S.hourlycarboused / self.pop_slab # original
        S.hourlycarboused = 0.
        # dividing it by PopSlab converts it to g/day/plant
        #ToDo: need to document this better, what is pcrs being used for.

        # SHOOTR->PCRS in 2dsoil is the actual rate of carbon supplied to roots in a soil slab, it is in g/day.
        # - dividing it by (SHOOTR->Rowsp*1)/10000 converts it to g/day/m^2
        # - further dividing it by weather->daylng converts it to g/hour/m^2
        # - then dividing it by plant density, converts it to g/hour/plant, which is the unit of the wthr.pcrs in maizesim. Yang. 10/27/06

        # Pass through nitrogen uptake (total mg per slab in the one hour) from 2DSOIL.
        w.TotalRootWeight = S.totalrootweight / self.pop_slab
        w.MaxRootDepth = S.maxrootdepth.item()
        # Available water is cm per profile - should be divided by PopSlab
        w.ThetaAvail = N.thetaavail / self.pop_slab

    def _update_weather_evapotranspiration(self, w, S):
        #FIXME LAI check might be not needed, as gas exchange module will deal with it anyways
        # ET_Supply is the actual amount of water that can be taken from the soil slab ( derived from AWUPS, g day-1 slab-1). To compare this variable with the et rate in maizesim it has to be converted into grams water per plant. To do this multiply by EOMULT to double slab width if plant is at the edge. Then multiply by 100/PopRow to get area inhabited by the plant. This provides a per plant estimate from area.
        if S.lai == 0:
            w.ET_supply = 0.
        else:
            # Note water uptake has been summed over the past hour so it is an hourly amount
            # into MAIZESIM Yang 8/15/06, dt 4/24/2011
            w.ET_supply = self.water_uptake / (S.eomult * S.poprow) * 100.
            #dt 4-24-2011 I replaced SHOOTR->AWUPS with WaterUptake. AWUPS is an instantaneous value.

            # for debugging
            ET_diff = w.ET_supply * 24 - S.et_demand

    def _update_plant_evapotranspiration(self, S):
        # ET_Supply is the actual amount of water that can be taken from the soil slab ( derived from AWUPS, g day-1 slab-1). To compare this variable with the et rate in maizesim it has to be converted into grams water per plant. To do this multiply by EOMULT to double slab width if plant is at the edge. Then multiply by 100/PopRow to get area inhabited by the plant. This provides a per plant estimate from area.
        if S.lai == 0:
            self.plant.water.supply = 0.
        else:
            # Note water uptake has been summed over the past hour so it is an hourly amount
            # into MAIZESIM Yang 8/15/06, dt 4/24/2011
            self.plant.water.supply = self.water_uptake / (S.eomult * S.poprow) * 100.
            #dt 4-24-2011 I replaced SHOOTR->AWUPS with WaterUptake. AWUPS is an instantaneous value.

            # for debugging
            ET_diff = self.plant.water.supply * 24 - S.et_demand

    def _update_weather_soil(self, w, G, N):
        # First find top of grid.
        y = G.y[:G.numnp]
        max_y = y.max()
        lower_boundary = max_y - 5.

        # Now find average temperature in layer between surface and lower boundary.
        ylb = (y >= lower_boundary)
        soil_t = N.tmpr[ylb].sum()
        count = ylb.sum()

        w.soilT = soil_t / (count*1.)

    ##############
    # Controller #
    ##############

    # The model code to simulate growth ect begins here when the plant object is called.
    def _run_controller(self):
        #TODO add some error catching code here
        # ier = self.controller.getErrStatus()
        # if ier != 0:
        #     raise Exception("controller is in error: " + ier)

        self.controller.run(
            Weather.from_2DSOIL(T, W),
            Soil.from_2DSOIL(T, W, S, N, G),
        )
        # Pass both weather and leaf water potential into the "run" function
        # of the controller pSC YY
        # need to get rid of other run module (with only wthr) done?
        # need to add PredawnLWP to wthr structure

    ###############
    # Development #
    ###############
    def _handle_not_germinated(self, w):
        # Assumes that germination takes place about halfway through the sowing date.
        if 0.49 <= w.time <= 0.51:
            print(w.jday)

    def _handle_emerged(self, w, S):
        # pass appropriate data to 2DSOIL file structures
        #dt 03/14/2011- I added a pool of carbo to hold leftover carbon from root growth, here it is implemented - see also plant
        # This holds any carbon not used for root growth in the previous time step
        plant = self.plant
        pool = plant.carbon.root_pool
        root = plant.carbon.root
        shoot = plant.carbon.shoot
        #ii = self.initinfo

        # this assures the pool is only used at night
        # minimizes complexity when pcrq has a value
        # since we have to add leftover carbo from pcrq to the shoot
        if pool > 0 and root < 0.00001:
            S.pcrl = (root + pool) * 24*self.pop_slab
            plant.carbon.reset_root_pool()
        else:
            S.pcrl = root * 24*self.pop_slab

        if self.plant.pheno.grain_filling:
            S.pcrq = (root + 0.75*shoot) * 24*self.pop_slab
        else:
            S.pcrq = (root + shoot) * 24*self.pop_slab

        #DT 09/19/14 under strong water stress mid season too much carbon is allocated to the roots,
        #we try to limit it here.
        #SHOOTR->PCRQ=SHOOTR->PCRL; //for debugging now remove later
        #dt 03/2011 added these two for debugging now
        #need to calculate mass balcance of carbo sent to root
        #can drop them later
        w.pcrl = S.pcrl / self.pop_slab/24.
        w.pcrq = S.pcrq / self.pop_slab/24.

        S.lcai = plant.area.green_leaf * self.initials.plant_density
        S.cover = 1 - np.exp(-0.79*S.lcai)
        S.shade = S.cover * S.rowsp
        S.height = min(S.shade, S.rowsp)
        S.et_demand = plant.water.supply * 24 # pass ET demand from shoot to root. Yang
        # In GasExchange, the unit of ET is mmol m-2(leaf) sec-1
        # need to convert to grams plant-1
        # Here, multiplying ET by 0.018 and 3600*24 converts it to g m-2(ground) day-1
        # dividing it by plantdensity converts it to g plant-1 day-1
        S.lai = S.lcai

        shoot_weight_per_m2 = plant.mass.shoot * self.initials.plant_density # Calculate total shoot mass per meter aquared YY
        mass_increase = shoot_weight_per_m2 - self.old_shoot_weight_per_m2 # Calculated increase in above-ground biomass per m2 YY

        # The biomass returned by getPlant()->get_shootMass() is the weight of each single plant (g/plant),
        # to convert it into (g/m-2), it has to be
        # multiplied by pSC->getIniInfo().plantDensity

        # Perhaps it's good idea to merge all plant N business into plant.cpp where C allocation is taken care of.
        # For now I'm not modifying any codes here. TODO - dt 03/15/2011 Much of this nitrogen stuff can be cleaned up as some of
        # Yang's original code is no longer used (i.e., U_N etc).

        #FIXME: move parameters out
        d = 0.075 # d: shape coefficient in the logistic function to simulate cumulative N uptake (Equation 9 in Lindquist et al. 2007)
        q_n = 0.032 # q_n: the maximum ratio of daily N uptake to measured daily growth rate (g N g-1) (Lindquist et al., 2007)
        a = 4.10 # maximum nitrogen concentration for C4 species: 4.1% (Lindquist et al. 2007)
        b = 0.5 # shape coefficient in calculation of nitrogen concentration in relation to up-ground biomass (equation 4 Lindquist et al, 2007) YY

        nitrogen_ratio = a / 100.
        if shoot_weight_per_m2 >= 100:
            # Calcualte above ground potential N concentration
            #nitrogen_ratio *= np.sqrt(shoot_weight_per_m2)
            nitrogen_ratio *= pow(shoot_weight_per_m2, 1 - b)
        #FIXME no nitrogen handling in Plant
        #plant.set_NitrogenRatio(nitrogen_ratio / 10.)

        # U_N maximum observed N uptake rate (g N m-2 ground d-1) (Lindquist et al, 2007) YY
        # The unit of U_N is g N m-2 ground d-1
        U_N = 0.359 * d / 4.

        # U_M maximum uptake rate as limited by maximum N fraction per unit (Equation 2 in Lindquist et al., 2007)
        # unit of U_M is also g N m-2 ground d-1; however, the unit of massIncrease is g m-2/step (here one hour).
        # Here in the simulation, the default length of one step is an hour; so, we have to scale it up to one
        # day by multiplying it by 24
        U_M = q_n * mass_increase * 24

        # unit of U_P is also g N m-2 ground d-1; however, the unit of massIncrease is g/step. Here
        # in the simulation, the default length of one step is an hour; so, we have to scale it up to one
        # day by multiplying it by 24

        # U_P potential rate of N accumulation (g N m-2 ground d-1) (Lindquist et al. 2007)
        # if shoot weight<100 (g m-2) then U_P is calculated this way
        U_P = (a / 100.) * mass_increase * 24
        # otherwise, it is calculated like this (Equation 6, Lindquist et al., 2007) YY
        if shoot_weight_per_m2 >= 100:
            U_P *= 10 * (1 - b) * pow(shoot_weight_per_m2, -b)

        # U_D U uptake rate (g N m-2 d-1) as limited by the difference between potential and actual amount of N
        # in existing biomass, equation 3 in Lindquist et al. 2007)
        # the returned value from get_N() is in g N/plant. It has to be converted to g/m-2 ground
        # that's why the actual n content is mulitpled by pSC->getIniInfo().plantDensity/(100*100) YY
        U_D = a * 10 / 100. * pow(shoot_weight_per_m2, -b) - self.plant.nitrogen.pool * self.initials.plant_density

        # set up account of N here
        # first hourly
        # Actual and needed N uptake in the last hour per plant per day

        # houly rate per day
        hourly_actual_n_from_soil = (self.nitrogen_uptake - self.nitrogen_uptake_old) / self.pop_slab
        #FIXME no nitrogen handling in Plant
        #plant.set_HourlyNitrogenSoilUptake(hourly_actual_n_from_soil)

        # Determine the nitrogen demand (equation 1 Lindquist et al. 2007) in grams plant-1
        hourly_nitrogen_demand = max(U_P, 0) / self.initials.plant_density / 24.
        #FIXME no nitrogen handling in Plant
        #plant.set_HourlyNitrogenDemand(hourly_nitrogen_demand)

        # now do cumulative amounts
        self.cumulative_nitrogen_demand += hourly_nitrogen_demand # grams plant-1 day-1
        #FIXME no nitrogen handling in Plant
        #plant.set_CumulativeNitrogenDemand(self.cumulative_nitrogen_demand)
        #plant.set_CumulativeNitrogenSoilUptake(self.nitrogen_uptake / self.pop_slab)

        # Pass the nitrogen demand into 2dsoil YY
        # Units are ug slab-1
        #oldNDemand = shootr.nitrodemand / popSlab / 1e6 / 24.
        S.nitrodemand = hourly_nitrogen_demand * self.pop_slab * 1e6 * 24

        # Save the value of the above_ground biomass of this time-step
        self.old_shoot_weight_per_m2 = shoot_weight_per_m2

        # save the cumulative N uptake from this time step
        self.nitrogen_uptake_old = self.nitrogen_uptake

        #HACK: moved to _update_nitrogen_uptake_error()
        #S.ndemanderror = self.nitrogen_demand_error
        #S.cumulativendemanderror = self.cumulative_nitrogen_demand_error

    def _handle_dead_or_not(self, M, T):
        if self.plant.pheno.dead:
            print("Completing crop simulation...")

            # tell 2dsoil that crops harvested
            M.nshoot = 0

            # set the next time so the model
            T.tnext[self.mod_num - 1] = 1e12

            # if matured points to nothing
            self.controller = None

            T.runflag = 0
        else:
            T.tnext[self.mod_num - 1] = T.time + self.period
            self.water_uptake = 0.
            #self.nitrogen_uptake = 0.

    ################
    # Post-process #
    ################

    def _postprocess(self):
        soil.carbon_partitioning()
        soil.rootgrow()
        soil.wateruptake()
        soil.soluteuptake()
        soil.gasuptake()

        soil.watermover()
        soil.solutemover()
        soil.heatmover()
        soil.gasmover()
        soil.soilnitrogen()
        soil.macrochem()
        #soil.massbl()

        if T.outputsoilyes > 0:
            soil.output()

    def run(self):
        self._init()

        #TODO: implement own synchronizer to properly exit the loop
        while True:
            self._preprocess()
            self._process()
            self._postprocess()

if __name__ == '__main__':
    driver = Driver()
    driver.run()
