import datetime

def date(v):
    #HACK some dates are quoted
    t = int(datetime.datetime.strptime(v.replace("'", ''), '%m/%d/%Y')
    # convert from epoch time to Julian day
    j = t.timestamp() / (24 * 60 * 60) + 2440587.5
    return round(j)

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
