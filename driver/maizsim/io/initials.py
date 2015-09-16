from .base import LegacyFile, date

class Initials(LegacyFile):
    @property
    def specs(self):
        return [
            [],
            [],
            [
                # PopRow: plant population per meter of row (m-1)
                ('population_per_row', float),

                # RowSP: row spacing (cm)
                ('row_spacing', float),

                # PopArea: plant density (m-2?)
                ('plant_density', float),

                # RowAng: row orientation measured eastward from North (degrees)
                ('row_angle', float),

                # xBStem: horizontal coordinate of the stem base (cm)
                ('stem_base_x', float),

                # yBStem: vertical coordinate of the stem base (cm)
                ('stem_base_y', float),

                # CEC: canopy extinction coefficient (0~1)
                ('canopy_extinction_coeff', float),

                # EOMult: multiplier depending on the plant position (0~1)
                # 1 if the plant is not on the border of the soil slab
                # 0 if it is on the border
                ('plant_unborderness', float),

                # Co2?
                #('co2', float),
            ],
            [],
            [
                # LATUDE: latitude (degrees)
                ('latitude', float),

                # Longitude: longitude (degrees)
                ('longitude', float),

                # Altitude: altitude (m)
                ('altitude', float),
            ],
            [],
            [
                # AutoIrrigate(F?)
                ('auto_irrigating', int),
            ],
            [],
            [
                # beginDay: simulation begin date (jday)
                ('begin_day', date),

                # sowingDay: sowing date (jday)
                ('sowing_day', date),

                # emergenceDay: sowing date (jday) -- old format?
                #('emergence_day', date),

                # endDay: simulation end date (jday)
                ('end_day', date),

                # TimeStep: timestep (min)
                ('timestep', int),
            ],
            [],
            [],
            [
                # OutputSoilNo (0/1)
                ('output_soil_no', int),

                # OutputSoilYes (0/1)
                ('output_soil_yes', int),
            ]
        ]
