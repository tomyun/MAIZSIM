from .base import LegacyFile, boolean

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
                ('day_length_sensitive', boolean),
                ('maximum_leaf_tip_appearance_rate', float), # Rmax_LTAR
                ('maximum_leaf_initiation_rate', float), # Rmax_LIR
                ('phyllochrons_to_silk', float),
            ]
        ]

        #TODO: support additional sections (2DSOIL exclusive): [SoilRoot], [SoilNitrogen]

    def __str__(self):
        return """\
Description: {description}
Cultivar: {cultivar}
GDD rating: {gdd_rating}
Generic leaf number: {generic_leaf_number}
Day length sensitive: {day_length_sensitive}
Maximum leaf tip appearance rate (Rmax_LTAR): {maximum_leaf_tip_appearance_rate}
Maximum leaf initiation rate (Rmax_LIR): {maximum_leaf_initiation_rate}
Phyllochrons to silk: {phyllochrons_to_silk}
""".format(**self.__dict__)
