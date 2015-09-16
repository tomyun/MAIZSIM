from .base import LegacyFile

class Variety(LegacyFile):
    @property
    def specs(self):
        return [
            [('description', str)],
            [('cultivar', str)],
            [],
            [],
            [
                ('gdd_rating', float),
                ('generic_leaf_number', int),
                ('day_length_sensitivity', int),
                ('Rmax_LTAR', float),
                ('Rmax_LTIR', float),
                ('phyllochrons_to_silk', float),
            ]
        ]
