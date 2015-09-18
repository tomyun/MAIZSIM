from .base import LegacyFile, date

class Time(LegacyFile):
    @property
    def specs(self):
        return [
            [],
            [],
            [
                # Time: initial time
                ('start_date', date),

                # Step: initial time increment dt (T)
                ('step', float),

                # dtMin: minimal acceptable time step (T)
                ('step_min', float),

                # dtMul1
                # If the number of iterations required at a particular time step is less than or equal to 3,
                # then dt for the next time step is multiplied by a dimensionless number dtMul1 >= 1.0
                ('dt_mul1', float),

                # dtMul2
                # If the number of iterations required is greater than or equal to 7,
                # then dt for the next time step is multiplied by dtMul2 < 7
                ('dt_mul2', float),

                # tFin: time to stop calculations
                ('end_date', date),
            ],
            [],
            [
                # output variables, 1 if true
                ('daily_output', int),
                ('hourly_output', int),
            ],
            [],
            [
                # weather data, 1 if true
                ('daily_weather', int),
                ('hourly_weather', int),
            ],
        ]
