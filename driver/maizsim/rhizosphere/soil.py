from ..timer import Timer

import copy

class Soil:
    def __init__(self):
        self.time = None
        self.T_soil = None
        self.ET_supply = None
        self.WP_leaf = None
        self.WP_leaf_predawn = None
        self.actual_root_carbon_supply_rate = None # PCRS (g day-1)
        self.min_root_carbon_supply_rate = None # PCRL (g day-1)
        self.max_root_carbon_supply_rate = None # PCRQ (g day-1)
        self.total_root_weight = None
        self.max_root_depth = None
        self.water = None # ThetaAvail (cm?)

    @classmethod
    def from_2DSOIL(cls, T, W, S, N, G):
        self.time = Timer.datetime_from_julian_day(T.time)
        self.update_water_potential(S)
        self.update_root(S, N)
        #FIXME handle ET_supply driectly in the driver
        #self.update_evapotranspiration(S)
        self.update_temperature(G, N)

    def update_water_potential(self, S):
        # since LeafWP in 2dsoil is in bar but in maizesim is in MPa,
        # so, have to divide it by 10 to convert it into MPa before passing the value to Maizesim 1 bar=10kPa
        self.WP_leaf = S.psil_ / 10. # and leaf water potential information into MAIZESIM Yang 8/15/06

        #FIXME make sure WP_leaf and WP_leaf_predawn have same unit
        # If time is 5 am, then pass the leaf water potential (the predawn leaf water potential) from SHOOTR to the wthr object. YY
        if self.time.hour == 5:
            # Here LeafWP is in bar.
            # Since the LWPeffect in leaf.cpp uses leaf water potential in bar,
            # so here PredawnLWP is in bar, instead of being scaled to MPa. YY
            self.WP_leaf_predawn = S.psil_.item() # SHOOTR->LeafWP

    def update_root(self, S, N):
        #plant_density = S.poprow * 100 / S.rowsp
        #return S.poprow / 100 * S.rowsp * S.eomult
        pop_slab = S.poprow / 100 * S.eomult

        self.min_root_carbon_supply_rate = S.pcrl / pop_slab / 24
        self.max_root_carbon_supply_rate = S.pcrq / pop_slab / 24

        # Pass actual carbohydrate amount used in 2dsoil back to the plant.
        #ToDo - make pcrs a new variable (ActualRootCarboUsed) and make it a member of plant.
        #dt here I changed this temporarily for debugging
        # don't need to divide by 24 since the value has been integrated over an hour
        self.actual_root_carbon_supply_rate = S.hourlycarboused / pop_slab # original
        #FIXME side-effect!
        S.hourlycarboused = 0.
        # dividing it by PopSlab converts it to g/day/plant
        #ToDo: need to document this better, what is pcrs being used for.

        # SHOOTR->PCRS in 2dsoil is the actual rate of carbon supplied to roots in a soil slab, it is in g/day.
        # - dividing it by (SHOOTR->Rowsp*1)/10000 converts it to g/day/m^2
        # - further dividing it by weather->daylng converts it to g/hour/m^2
        # - then dividing it by plant density, converts it to g/hour/plant, which is the unit of the wthr.pcrs in maizesim. Yang. 10/27/06

        # Pass through nitrogen uptake (total mg per slab in the one hour) from 2DSOIL.
        self.total_root_weight = S.totalrootweight / pop_slab
        self.max_root_depth = S.maxrootdepth.item()

        # Available water is cm per profile - should be divided by PopSlab
        self.water = N.thetaavail / pop_slab

    def update_temperature(self, G, N):
        # First find top of grid.
        y = G.y[:G.numnp]
        max_y = y.max()
        lower_boundary = max_y - 5

        # Now find average temperature in layer between surface and lower boundary.
        ylb = (y >= lower_boundary)
        soil_t = N.tmpr[ylb].sum()
        count = ylb.sum()
        self.T_soil = soil_t / count

    def copy(self):
        return copy.copy(self)
