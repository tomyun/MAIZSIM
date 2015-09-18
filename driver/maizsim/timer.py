import datetime

JULIAN_EPOCH_USDA = 2415078.5 # 1990-03-01
JULIAN_EPOCH_UNIX = 2440587.5 # 1970-01-01

class Timer:
    def __init__(self, time, step):
        self.time = time
        self.step = step

    @classmethod
    def from_datetime(cls, time, step):
        return cls(time, step)

    @classmethod
    def from_julian_day(cls, jday, step):
        time = cls.datetime_from_julian_day(jday)
        return cls.from_datetime(time, step)

    @staticmethod
    def datetime_from_julian_day(jday):
        d = jday + (JULIAN_EPOCH_USDA - JULIAN_EPOCH_UNIX)
        return datetime.datetime.utcfromtimestamp(d * (24 * 60 * 60))

    @staticmethod
    def julian_day_from_datetime(time):
        j = time.timestamp() / (24 * 60 * 60) - (JULIAN_EPOCH_USDA - JULIAN_EPOCH_UNIX)
        return round(j)

    def tick(self):
        self.time += datetime.timedelta(hours=self.step)