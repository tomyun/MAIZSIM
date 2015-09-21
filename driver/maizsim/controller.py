from .timer import Timer
from .physiology import Plant

import pandas as pd

class Controller:
    def __init__(self, initials, variety, time, crop_filename, leaf_filename):
        self.initials = initials
        self.variety = variety
        self.time = time
        self.crop_filename = crop_filename
        self.leaf_filename = leaf_filename
        self.setup()
        self.setup_crop_output()
        self.setup_leaf_output()

    def setup(self):
        self.show_title()
        self.show_variety()
        # Timer class gets stepsize in hours
        self.timer = Timer.from_julian_day(self.initials.begin_day, self.initials.timestep / 60)
        #TODO load up atmosphere.Weather here?
        #dt modified this after modifying code to use julian day
        #self.weather =
        self.plant = Plant(self.initials, self.variety)

    def show_title(self):
        print("""\
***********************************************************
*         PyMAIZSIM: A Simulation Model for Corn          *
*                     VERSION  1.1.00 2015                *
*    UNIVERSITY OF WASHINGTON, PLANT ECOPHYSIOLOGY LAB    *
*   USDA-ARS, CROP SYSTEMS AND GLOBAL CHANGE LABORATORY   *
***********************************************************
""")

    def show_variety(self):
        print(self.variety)

        sowing_day = Timer.datetime_from_julian_day(self.initials.sowing_day)
        print("Sowing day: {}".format(sowing_day))

    def run(self, weather, soil):
        #FIXME no need to check range here, since the current weather is updated every iteration?
        #if (weather[iCur].jday >= initInfo.sowingDay && weather[iCur].jday <= lastDayOfSim)
        self.plant.update(weather, soil)

        #dt added ability to output daily based on 2dsoil output
        # Always hourly for now - have to add code to average values
        if self.time.daily_output == 1 and weather.time.hour == 6:
            self.update_crop_output(weather, soil)
            if self.plant.pheno.germinated:
                self.update_leaf_output(weather, soil)
        elif self.time.hourly_output == 1:
            self.update_crop_output(weather, soil)
            if self.plant.pheno.germinated:
                self.update_leaf_output(weather, soil)

        # update timer step forward
        self.timer.tick()

    def setup_crop_output(self):
        columns = [
            'date',
            'jday',
            'time',
            'Leaves',
            'Dropped',
            'LA/pl',
            'LA_dead',
            'LAI',
            'RH',
            'LeafWP',
            'PFD',
            'SolRad',
            'SoilT',
            'Tair',
            'Tcan',
            'ETdmd',
            'ETsply',
            'Pn',
            'Pg',
            'Respir',
            'av_gs',
            'VPD',
            'Nitr',
            'N_Dem',
            'NUpt',
            'LeafN',
            'PCRL',
            'totalDM',
            'shootDM',
            'earDM',
            'GrleafDM',
            'DrpLfDM',
            'stemDM',
            'rootDM',
            'SoilRt',
            'MxRtDep',
            'AvailW',
            'solubleC',
            'Note',
        ]
        self.crop_output = pd.DataFrame(columns=columns)

    def update_crop_output(self, weather, soil):
        self.crop_output.loc[len(self.crop_output)] = [
            weather.time,
            Timer.julian_day_from_datetime(weather.time),
            weather.time.hour,
            self.plant.pheno.leaves_appeared,
            self.plant.count.total_dropped_leaves,
            self.plant.area.green_leaf,
            self.plant.area.senescent_leaf,
            self.plant.area.leaf_area_index,
            weather.RH,
            #None, # self.plant.carbon_ratio,
            soil.WP_leaf, # print out leaf water potential Yang 8/22/06
            weather.PFD,
            None, # weather.sol_rad,
            soil.T_soil,
            weather.T_air,
            self.plant.pheno.temperature,
            None, # self.plant.water.get_ET_old()
            self.plant.water.supply, # in both cases transpiration is grams per plant per hour
            self.plant.photosynthesis.net, # g Carbo per plant per hour
            self.plant.photosynthesis.gross,
            self.plant.carbon.maintenance_respiration, #dt 03/2011 added to better calc mass balance
            self.plant.photosynthesis.conductance if self.plant.pheno.emerged else 0, # return average stomatal conductance Yang 10/31/06
            self.plant.photosynthesis.vapor_pressure_deficit,
            self.plant.nitrogen.pool,
            None, # self.plant.nitrogen.CumulativeNitrogenDemand()
            None, # self.plant.nitrogen.CumulativeNitrogenSoilUptake()
            self.plant.nitrogen.leaf, # return mass of N in leaves YY
            None, # self.plant.carbon.ActualCarboIncrement()
            self.plant.mass.total,
            self.plant.mass.shoot, # masses are grams per plant
            self.plant.mass.ear,
            self.plant.mass.leaf,
            self.plant.mass.dropped_leaf,
            self.plant.mass.stem,
            self.plant.mass.root,
            soil.total_root_weight,
            soil.max_root_depth,
            soil.water,
            self.plant.carbon.reserve,
            self.plant.pheno.current_stage,
            #weather.day_length, #FIXME no column defined
            #self.plant.pheno.gdd_after_emergence,
        ]

    def export_crop_output(self):
        self.crop_ouptut.to_csv(self.crop_filename)

    def setup_leaf_output(self):
        columns = [
            'date',
            'jday',
            'time',
            'Lvs_Init',
            'Lvs_Apr',
            'Leaf_#',
            'area',
            'mass',
            'Sen_Area',
            'Pntl_Area',
            'Longev',
            'CarbRat',
            'SLA',
            'dropped',
            'state',
            'GDD Sum',
        ]
        self.leaf_output = pd.DataFrame(columns=columns)

    def update_leaf_output(self, weather, soil):
        def row(leaf):
            return [
                weather.time,
                Timer.julian_day_from_datetime(weather.time),
                weather.time.hour,
                self.plant.pheno.leaves_initiated,
                self.plant.pheno.leaves_appeared,
                leaf.rank,
                leaf.green_area,
                leaf.mass,
                leaf.senescent_area,
                leaf.potential_area,
                leaf.longevity,
                leaf.nitrogen,
                leaf.specific_leaf_area,
                leaf.dropped,
                leaf.growing,
                self.plant.pheno.gdd_after_emergence,
            ]
        for nu in self.plant.nodal_units:
            self.leaf_output.loc[len(self.leaf_output)] = row(nu.leaf)

    def export_leaf_output(self):
        self.leaf_ouptut.to_csv(self.leaf_filename)
