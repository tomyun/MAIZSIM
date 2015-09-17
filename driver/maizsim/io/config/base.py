import datetime
from ...timer import Timer

def date(v):
    #HACK some dates are quoted
    t = datetime.datetime.strptime(v.replace("'", ''), '%m/%d/%Y')
    return Timer.julian_day_from_datetime(t)

class LegacyFile:
    def __init__(self, filename):
        self.load(filename)

    @property
    def specs(self):
        return [[]]

    def load(self, filename):
        with open(filename) as f:
            self._parse(f)

    def _parse(self, f):
        for spec in self.specs:
            l = f.readline().strip()
            if len(spec) > 1:
                values = l.split()
            else:
                values = [l]
            [setattr(self, n, c(v)) for (n, c), v in zip(spec, values)]
