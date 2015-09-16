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
                ('maximum_leaf_tip_appearance_rate', float), # Rmax_LTAR
                ('maximum_leaf_initiation_rate', float), # Rmax_LIR
                ('phyllochrons_to_silk', float),
            ]
        ]

        #TODO: support additional sections (2DSOIL exclusive): [SoilRoot], [SoilNitrogen]
