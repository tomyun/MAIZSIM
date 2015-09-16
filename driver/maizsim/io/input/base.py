import datetime

def date(v):
    #HACK some dates are quoted
    return int(datetime.datetime.strptime(v.replace("'", ''), '%m/%d/%Y').strftime('%j'))

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
